use std::{
    fs,
    io::{self, Write},
    net::{Ipv4Addr, SocketAddr, SocketAddrV4},
    path::{Path, PathBuf},
    sync::{
        Arc, Condvar, Mutex, RwLock,
        atomic::{AtomicBool, Ordering},
    },
    thread::{self, JoinHandle},
    time::Duration,
};

use eyre::{Result, WrapErr, eyre};
use percent_encoding::percent_decode_str;
use tempfile::TempDir;
use tiny_http::{Header, Method, Request, Response, Server, StatusCode};

const RELOAD_SCRIPT: &str = r#"<script data-alphal00p-live-reload>(()=>{const events=new EventSource("/__events");events.addEventListener("reload",()=>location.reload())})();</script>"#;
const SSE_HEARTBEAT_INTERVAL: Duration = Duration::from_secs(15);

struct PublishedGeneration {
    root: PathBuf,
    _directory: TempDir,
}

#[derive(Default)]
struct ReloadState {
    generation: u64,
    stopped: bool,
}

#[derive(Default)]
struct ReloadHub {
    state: Mutex<ReloadState>,
    changed: Condvar,
}

#[derive(Debug, Eq, PartialEq)]
enum ReloadWait {
    Published(u64),
    Heartbeat,
    Stopped,
}

impl ReloadHub {
    fn generation(&self) -> u64 {
        self.state.lock().expect("reload state poisoned").generation
    }

    fn publish(&self) {
        let mut state = self.state.lock().expect("reload state poisoned");
        state.generation += 1;
        self.changed.notify_all();
    }

    fn wait_after(&self, generation: u64, timeout: Duration) -> ReloadWait {
        let state = self.state.lock().expect("reload state poisoned");
        let (state, wait) = self
            .changed
            .wait_timeout_while(state, timeout, |state| {
                !state.stopped && state.generation <= generation
            })
            .expect("reload state poisoned");
        if state.stopped {
            ReloadWait::Stopped
        } else if state.generation > generation {
            ReloadWait::Published(state.generation)
        } else {
            debug_assert!(wait.timed_out());
            ReloadWait::Heartbeat
        }
    }

    fn stop(&self) {
        let mut state = self.state.lock().expect("reload state poisoned");
        state.stopped = true;
        self.changed.notify_all();
    }
}

struct ServerState {
    publication: RwLock<Option<Arc<PublishedGeneration>>>,
    status: RwLock<String>,
    reload: Arc<ReloadHub>,
}

impl Default for ServerState {
    fn default() -> Self {
        Self {
            publication: RwLock::new(None),
            status: RwLock::new("Waiting for the first successful build.".to_owned()),
            reload: Arc::new(ReloadHub::default()),
        }
    }
}

pub(crate) struct LiveServer {
    address: SocketAddr,
    state: Arc<ServerState>,
    server: Arc<Server>,
    running: Arc<AtomicBool>,
    thread: Option<JoinHandle<()>>,
}

impl LiveServer {
    pub(crate) fn start(port: u16) -> Result<Self> {
        let server = Arc::new(
            Server::http(SocketAddrV4::new(Ipv4Addr::LOCALHOST, port))
                .map_err(|error| eyre!("failed to bind documentation server: {error}"))?,
        );
        let address = server
            .server_addr()
            .to_ip()
            .ok_or_else(|| eyre!("documentation server did not bind a TCP address"))?;
        let state = Arc::new(ServerState::default());
        let running = Arc::new(AtomicBool::new(true));
        let thread = {
            let server = Arc::clone(&server);
            let state = Arc::clone(&state);
            let running = Arc::clone(&running);
            thread::spawn(move || {
                while running.load(Ordering::Acquire) {
                    match server.recv_timeout(Duration::from_millis(100)) {
                        Ok(Some(request)) => {
                            let state = Arc::clone(&state);
                            thread::spawn(move || {
                                if let Err(error) = handle_request(request, &state) {
                                    eprintln!("documentation server request failed: {error:#}");
                                }
                            });
                        }
                        Ok(None) => {}
                        Err(error) if running.load(Ordering::Acquire) => {
                            eprintln!("documentation server failed: {error}");
                        }
                        Err(_) => break,
                    }
                }
            })
        };
        Ok(Self {
            address,
            state,
            server,
            running,
            thread: Some(thread),
        })
    }

    pub(crate) fn address(&self) -> SocketAddr {
        self.address
    }

    pub(crate) fn set_status(&self, status: impl Into<String>) {
        *self.state.status.write().expect("server status poisoned") = status.into();
    }

    pub(crate) fn publish(&self, directory: TempDir, root: PathBuf) -> Result<()> {
        if !root.join("index.html").is_file() {
            return Err(eyre!(
                "documentation generation has no repository index at {}",
                root.display()
            ));
        }
        let generation = Arc::new(PublishedGeneration {
            root,
            _directory: directory,
        });
        *self
            .state
            .publication
            .write()
            .expect("publication state poisoned") = Some(generation);
        self.set_status("Serving the last successful build.");
        self.state.reload.publish();
        Ok(())
    }
}

impl Drop for LiveServer {
    fn drop(&mut self) {
        self.running.store(false, Ordering::Release);
        self.state.reload.stop();
        self.server.unblock();
        if let Some(thread) = self.thread.take() {
            let _ = thread.join();
        }
    }
}

fn handle_request(request: Request, state: &Arc<ServerState>) -> Result<()> {
    if request.url().split('?').next() == Some("/__events") {
        return serve_events(request, Arc::clone(&state.reload));
    }
    if !matches!(request.method(), Method::Get | Method::Head) {
        return respond(
            request,
            StatusCode(405),
            "text/plain; charset=utf-8",
            b"method not allowed".to_vec(),
            None,
        );
    }

    let publication = state
        .publication
        .read()
        .expect("publication state poisoned")
        .clone();
    let Some(publication) = publication else {
        let status = state.status.read().expect("server status poisoned");
        let body = format!(
            "<!doctype html><meta charset=\"utf-8\"><meta http-equiv=\"refresh\" content=\"2\"><title>Documentation build</title><style>body{{max-width:45rem;margin:10vh auto;padding:2rem;font:18px/1.5 system-ui;color:#172033}}</style><h1>Documentation build</h1><p>{}</p>{RELOAD_SCRIPT}",
            escape_html(&status),
        );
        return respond(
            request,
            StatusCode(200),
            "text/html; charset=utf-8",
            body.into_bytes(),
            None,
        );
    };

    let url_path = request.url().split('?').next().unwrap_or("/").to_owned();
    let relative = match safe_relative_url_path(&url_path) {
        Ok(path) => path,
        Err(_) => {
            return respond(
                request,
                StatusCode(400),
                "text/plain; charset=utf-8",
                b"bad request".to_vec(),
                None,
            );
        }
    };
    let mut path = publication.root.join(&relative);
    if path.is_dir() {
        if !url_path.ends_with('/') {
            let location = format!("{url_path}/");
            return respond(
                request,
                StatusCode(308),
                "text/plain; charset=utf-8",
                Vec::new(),
                Some(("Location", &location)),
            );
        }
        path.push("index.html");
    }
    if !path.is_file() {
        return respond(
            request,
            StatusCode(404),
            "text/plain; charset=utf-8",
            b"not found".to_vec(),
            None,
        );
    }

    let content_type = content_type(&path);
    let mut body =
        fs::read(&path).wrap_err_with(|| format!("failed to read {}", path.display()))?;
    if content_type.starts_with("text/html") {
        let html = String::from_utf8(body).wrap_err("generated HTML is not UTF-8")?;
        body = inject_live_reload(&html).into_bytes();
    }
    respond(request, StatusCode(200), content_type, body, None)
}

fn serve_events(request: Request, hub: Arc<ReloadHub>) -> Result<()> {
    if request.method() != &Method::Get {
        return respond(
            request,
            StatusCode(405),
            "text/plain; charset=utf-8",
            b"method not allowed".to_vec(),
            None,
        );
    }
    // tiny_http buffers its chunked response writer until the reader reaches EOF. SSE never
    // reaches EOF while it is useful, so take ownership of the socket writer and flush each
    // event ourselves.
    let mut seen = hub.generation();
    let mut writer = request.into_writer();
    write_event_bytes(
        &mut writer,
        b"HTTP/1.1 200 OK\r\nContent-Type: text/event-stream\r\nCache-Control: no-cache\r\n\r\n: connected\n\n",
    )
    .wrap_err("failed to establish live-reload event stream")?;
    loop {
        let event = match hub.wait_after(seen, SSE_HEARTBEAT_INTERVAL) {
            ReloadWait::Published(generation) => {
                seen = generation;
                b"event: reload\ndata:\n\n".as_slice()
            }
            ReloadWait::Heartbeat => b": heartbeat\n\n".as_slice(),
            ReloadWait::Stopped => break,
        };
        if let Err(error) = write_event_bytes(&mut writer, event) {
            return match error.kind() {
                io::ErrorKind::BrokenPipe
                | io::ErrorKind::ConnectionAborted
                | io::ErrorKind::ConnectionReset => Ok(()),
                _ => Err(error).wrap_err("failed to send live-reload event"),
            };
        }
    }
    Ok(())
}

fn write_event_bytes(writer: &mut dyn Write, bytes: &[u8]) -> io::Result<()> {
    writer.write_all(bytes)?;
    writer.flush()
}

fn respond(
    request: Request,
    status: StatusCode,
    content_type: &str,
    body: Vec<u8>,
    extra_header: Option<(&str, &str)>,
) -> Result<()> {
    let mut response = Response::from_data(body)
        .with_status_code(status)
        .with_header(header("Content-Type", content_type))
        .with_header(header("Cache-Control", "no-store"));
    if let Some((name, value)) = extra_header {
        response.add_header(header(name, value));
    }
    request
        .respond(response)
        .wrap_err("failed to send documentation response")
}

fn header(name: &str, value: &str) -> Header {
    Header::from_bytes(name.as_bytes(), value.as_bytes()).expect("valid HTTP header")
}

fn safe_relative_url_path(url: &str) -> Result<PathBuf> {
    let decoded = percent_decode_str(url)
        .decode_utf8()
        .wrap_err("request path is not valid UTF-8")?;
    let mut path = PathBuf::new();
    for segment in decoded.split('/') {
        if segment.is_empty() {
            continue;
        }
        if segment == "." || segment == ".." || segment.contains(['\\', '\0']) {
            return Err(eyre!("unsafe documentation request path"));
        }
        path.push(segment);
    }
    Ok(path)
}

fn content_type(path: &Path) -> &'static str {
    match path.extension().and_then(|extension| extension.to_str()) {
        Some("html") => "text/html; charset=utf-8",
        Some("css") => "text/css; charset=utf-8",
        Some("js") => "text/javascript; charset=utf-8",
        Some("json") => "application/json",
        Some("pdf") => "application/pdf",
        Some("svg") => "image/svg+xml",
        Some("png") => "image/png",
        Some("jpg" | "jpeg") => "image/jpeg",
        Some("woff") => "font/woff",
        Some("woff2") => "font/woff2",
        Some("ttf") => "font/ttf",
        Some("txt" | "pyi") => "text/plain; charset=utf-8",
        _ => "application/octet-stream",
    }
}

fn inject_live_reload(html: &str) -> String {
    if html.contains("data-alphal00p-live-reload") {
        return html.to_owned();
    }
    if let Some(body) = html.rfind("</body>") {
        let mut output = String::with_capacity(html.len() + RELOAD_SCRIPT.len());
        output.push_str(&html[..body]);
        output.push_str(RELOAD_SCRIPT);
        output.push_str(&html[body..]);
        output
    } else {
        format!("{html}{RELOAD_SCRIPT}")
    }
}

fn escape_html(value: &str) -> String {
    value
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

#[cfg(test)]
mod tests {
    use std::{
        io::{BufReader, Read, Write},
        net::{Shutdown, TcpStream},
        time::Duration,
    };

    use super::*;

    struct TestResponse {
        status: u16,
        headers: Vec<(String, String)>,
        body: Vec<u8>,
    }

    impl TestResponse {
        fn header(&self, name: &str) -> Option<&str> {
            self.headers
                .iter()
                .find(|(candidate, _)| candidate.eq_ignore_ascii_case(name))
                .map(|(_, value)| value.as_str())
        }
    }

    fn request_with_method(address: SocketAddr, method: &str, target: &str) -> TestResponse {
        let mut stream = TcpStream::connect(address).unwrap();
        stream
            .set_read_timeout(Some(Duration::from_secs(5)))
            .unwrap();
        write!(
            stream,
            "{method} {target} HTTP/1.1\r\nHost: localhost\r\nConnection: close\r\n\r\n"
        )
        .unwrap();
        stream.shutdown(Shutdown::Write).unwrap();

        let mut response = Vec::new();
        stream.read_to_end(&mut response).unwrap();
        let header_end = response
            .windows(4)
            .position(|window| window == b"\r\n\r\n")
            .expect("HTTP response has a header terminator");
        let head = String::from_utf8(response[..header_end].to_vec()).unwrap();
        let mut lines = head.lines();
        let status = lines
            .next()
            .unwrap()
            .split_whitespace()
            .nth(1)
            .unwrap()
            .parse()
            .unwrap();
        let headers = lines
            .map(|line| {
                let (name, value) = line.split_once(':').unwrap();
                (name.to_owned(), value.trim().to_owned())
            })
            .collect();
        TestResponse {
            status,
            headers,
            body: response[header_end + 4..].to_vec(),
        }
    }

    fn request(address: SocketAddr, target: &str) -> TestResponse {
        request_with_method(address, "GET", target)
    }

    fn generation(index: &str) -> (TempDir, PathBuf) {
        let directory = tempfile::tempdir().unwrap();
        let root = directory.path().join("site");
        fs::create_dir(&root).unwrap();
        fs::write(root.join("index.html"), index).unwrap();
        (directory, root)
    }

    fn read_until(reader: &mut BufReader<TcpStream>, needle: &[u8]) -> Vec<u8> {
        let mut received = Vec::new();
        let mut buffer = [0; 256];
        while !received
            .windows(needle.len())
            .any(|window| window == needle)
        {
            let count = reader.read(&mut buffer).unwrap();
            assert_ne!(count, 0, "event stream ended before the expected message");
            received.extend_from_slice(&buffer[..count]);
        }
        received
    }

    #[test]
    fn live_reload_is_injected_once() {
        let html = "<html><body>hello</body></html>";
        let once = inject_live_reload(html);
        let twice = inject_live_reload(&once);
        assert_eq!(once, twice);
        assert_eq!(once.matches("data-alphal00p-live-reload").count(), 1);
        assert!(once.contains("EventSource(\"/__events\")"));
    }

    #[test]
    fn request_paths_are_decoded_without_allowing_traversal() {
        assert_eq!(
            safe_relative_url_path("/products/Spenso%20Docs/").unwrap(),
            Path::new("products/Spenso Docs")
        );
        assert!(safe_relative_url_path("/../Cargo.toml").is_err());
        assert!(safe_relative_url_path("/%2e%2e/Cargo.toml").is_err());
        assert!(safe_relative_url_path("/docs%5csecret").is_err());
    }

    #[test]
    fn server_shows_the_build_status_before_the_first_publication() {
        let server = LiveServer::start(0).unwrap();
        server.set_status("Compiling <guides> & API references.");

        let response = request(server.address(), "/");
        assert_eq!(response.status, 200);
        assert_eq!(
            response.header("Content-Type"),
            Some("text/html; charset=utf-8")
        );
        let body = String::from_utf8(response.body).unwrap();
        assert!(body.contains("Documentation build"));
        assert!(body.contains("Compiling &lt;guides&gt; &amp; API references."));
        assert!(body.contains("http-equiv=\"refresh\""));
        assert!(body.contains("EventSource(\"/__events\")"));
    }

    #[test]
    fn publication_injects_html_but_serves_assets_byte_for_byte() {
        let server = LiveServer::start(0).unwrap();
        let (directory, root) = generation("<!doctype html><html><body><p>ready</p></body></html>");
        let asset = [0, 255, 1, 254, b'\r', b'\n', 42];
        fs::create_dir(root.join("assets")).unwrap();
        fs::write(root.join("assets/example.bin"), asset).unwrap();
        server.publish(directory, root).unwrap();

        let page = request(server.address(), "/");
        assert_eq!(page.status, 200);
        let html = String::from_utf8(page.body).unwrap();
        assert!(html.contains("<p>ready</p>"));
        assert_eq!(html.matches("data-alphal00p-live-reload").count(), 1);
        assert!(html.find(RELOAD_SCRIPT).unwrap() < html.find("</body>").unwrap());

        let binary = request(server.address(), "/assets/example.bin");
        assert_eq!(binary.status, 200);
        assert_eq!(
            binary.header("Content-Type"),
            Some("application/octet-stream")
        );
        assert_eq!(binary.body, asset);
    }

    #[test]
    fn head_reports_the_get_representation_length_without_a_body() {
        let server = LiveServer::start(0).unwrap();
        let (directory, root) = generation("<html><body>ready</body></html>");
        server.publish(directory, root).unwrap();

        let get = request(server.address(), "/");
        let head = request_with_method(server.address(), "HEAD", "/");
        assert_eq!(head.status, get.status);
        assert_eq!(head.header("Content-Type"), get.header("Content-Type"));
        assert_eq!(head.header("Content-Length"), get.header("Content-Length"));
        assert_eq!(
            head.header("Content-Length"),
            Some(get.body.len().to_string().as_str())
        );
        assert!(head.body.is_empty());
    }

    #[test]
    fn failed_status_keeps_the_last_good_generation() {
        let server = LiveServer::start(0).unwrap();
        let (first_directory, first_root) =
            generation("<html><body>first generation</body></html>");
        server.publish(first_directory, first_root).unwrap();

        server.set_status("The latest build failed.");
        let retained = request(server.address(), "/");
        assert_eq!(retained.status, 200);
        let retained = String::from_utf8(retained.body).unwrap();
        assert!(retained.contains("first generation"));
        assert!(!retained.contains("latest build failed"));

        let (second_directory, second_root) =
            generation("<html><body>second generation</body></html>");
        server.publish(second_directory, second_root).unwrap();
        let updated = String::from_utf8(request(server.address(), "/").body).unwrap();
        assert!(updated.contains("second generation"));
        assert!(!updated.contains("first generation"));
    }

    #[test]
    fn routes_redirect_directories_and_reject_missing_or_unsafe_paths() {
        let server = LiveServer::start(0).unwrap();
        let (directory, root) = generation("<html><body>home</body></html>");
        fs::create_dir(root.join("guide")).unwrap();
        fs::write(
            root.join("guide/index.html"),
            "<html><body>guide</body></html>",
        )
        .unwrap();
        server.publish(directory, root).unwrap();

        let redirect = request(server.address(), "/guide");
        assert_eq!(redirect.status, 308);
        assert_eq!(redirect.header("Location"), Some("/guide/"));

        let guide = request(server.address(), "/guide/");
        assert_eq!(guide.status, 200);
        assert!(String::from_utf8(guide.body).unwrap().contains("guide"));

        let missing = request(server.address(), "/does-not-exist");
        assert_eq!(missing.status, 404);
        assert_eq!(missing.body, b"not found");

        let traversal = request(server.address(), "/%2e%2e/Cargo.toml");
        assert_eq!(traversal.status, 400);
        assert_eq!(traversal.body, b"bad request");
    }

    #[test]
    fn publication_notifies_an_http_event_stream() {
        let server = LiveServer::start(0).unwrap();
        let mut stream = TcpStream::connect(server.address()).unwrap();
        stream
            .set_read_timeout(Some(Duration::from_secs(5)))
            .unwrap();
        write!(
            stream,
            "GET /__events HTTP/1.1\r\nHost: localhost\r\nConnection: close\r\n\r\n"
        )
        .unwrap();
        let mut reader = BufReader::new(stream);
        let connected = read_until(&mut reader, b": connected\n\n");
        assert!(
            connected
                .windows(12)
                .any(|window| window == b"HTTP/1.1 200")
        );
        let connected = String::from_utf8(connected).unwrap();
        assert!(connected.contains("Content-Type: text/event-stream"));
        assert!(!connected.to_ascii_lowercase().contains("transfer-encoding"));

        let (directory, root) = generation("<html><body>ready</body></html>");
        server.publish(directory, root).unwrap();
        read_until(&mut reader, b"event: reload\ndata:\n\n");
    }

    #[test]
    fn reload_hub_wakes_every_connected_client() {
        let hub = Arc::new(ReloadHub::default());
        let generation = hub.generation();
        let readers = (0..2)
            .map(|_| {
                let hub = Arc::clone(&hub);
                thread::spawn(move || hub.wait_after(generation, Duration::from_secs(1)))
            })
            .collect::<Vec<_>>();
        hub.publish();
        for reader in readers {
            assert_eq!(reader.join().unwrap(), ReloadWait::Published(1));
        }
    }

    #[test]
    fn reload_hub_heartbeats_while_idle_and_stops_cleanly() {
        let hub = ReloadHub::default();
        assert_eq!(
            hub.wait_after(hub.generation(), Duration::ZERO),
            ReloadWait::Heartbeat
        );
        hub.stop();
        assert_eq!(
            hub.wait_after(hub.generation(), Duration::ZERO),
            ReloadWait::Stopped
        );
    }
}
