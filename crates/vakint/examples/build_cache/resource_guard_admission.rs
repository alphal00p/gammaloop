//! Optional task-scoped admission before Native work; policy belongs to the reviewed helper.

use std::{
    ffi::OsString,
    fs::OpenOptions,
    io::{self, Read},
    path::{Path, PathBuf},
    process::{Command, Stdio},
};

use serde_json::Value;
use sha2::{Digest, Sha256};

const ENVIRONMENT: [&str; 5] = [
    "VAKINT_RESOURCE_GUARD_PYTHON",
    "VAKINT_RESOURCE_GUARD_HELPER",
    "VAKINT_RESOURCE_GUARD_HELPER_SHA256",
    "VAKINT_RESOURCE_GUARD_LAUNCH_SPEC",
    "VAKINT_RESOURCE_GUARD_LAUNCH_SPEC_SHA256",
];
const MAX_METADATA_BYTES: u64 = 1024 * 1024;

pub(super) fn admit_from_task_env() -> io::Result<Option<Value>> {
    let Some(request) = Request::from_values(ENVIRONMENT.map(std::env::var_os))? else {
        // No process, file or helper access for an ordinary unconfigured call.
        return Ok(None);
    };
    let worker = WorkerIdentity::current()?;
    request.admit(&worker).map(Some)
}

struct Request {
    python: PathBuf,
    helper: PathBuf,
    helper_sha256: String,
    launch_spec: PathBuf,
    launch_spec_sha256: String,
}

impl Request {
    fn from_values(values: [Option<OsString>; 5]) -> io::Result<Option<Self>> {
        if values.iter().all(Option::is_none) {
            return Ok(None);
        }
        let [
            Some(python),
            Some(helper),
            Some(helper_sha256),
            Some(launch_spec),
            Some(launch_spec_sha256),
        ] = values
        else {
            return Err(io::Error::other(
                "all five VAKINT_RESOURCE_GUARD_* admission variables are required",
            ));
        };
        let hash = |value: OsString| -> io::Result<String> {
            let value = value
                .into_string()
                .map_err(|_| io::Error::other("guard hashes must be UTF-8"))?;
            if !valid_sha256(&value) {
                return Err(io::Error::other("guard hashes must be lowercase SHA256"));
            }
            Ok(value)
        };
        let request = Self {
            python: python.into(),
            helper: helper.into(),
            helper_sha256: hash(helper_sha256)?,
            launch_spec: launch_spec.into(),
            launch_spec_sha256: hash(launch_spec_sha256)?,
        };
        if [&request.python, &request.helper, &request.launch_spec]
            .iter()
            .any(|path| !path.is_absolute())
        {
            return Err(io::Error::other(
                "guard interpreter, helper and launch specification must be absolute paths",
            ));
        }
        Ok(Some(request))
    }

    fn verify_inputs(&self) -> io::Result<()> {
        if !self.python.is_file()
            || fingerprint(&self.helper, Some(MAX_METADATA_BYTES))? != self.helper_sha256
            || fingerprint(&self.launch_spec, Some(MAX_METADATA_BYTES))? != self.launch_spec_sha256
        {
            return Err(io::Error::other(
                "guard interpreter/helper/launch specification binding differs",
            ));
        }
        Ok(())
    }

    fn admit(&self, worker: &WorkerIdentity) -> io::Result<Value> {
        self.verify_inputs()?;
        // No shell and no policy reimplementation. The helper authenticates its
        // direct Rust parent and the authoritative supervisor/unit/ticket, then
        // waits for the existing guard's exact PIDFD acknowledgement and latches.
        let output = Command::new(&self.python)
            .arg("-B")
            .arg(&self.helper)
            .arg("--direct-parent")
            .arg("--launch-spec")
            .arg(&self.launch_spec)
            .arg("--launch-spec-sha256")
            .arg(&self.launch_spec_sha256)
            .stdin(Stdio::null())
            .stderr(Stdio::inherit())
            .output()?;
        // The helper's final action authenticates a fresh PIDFD acknowledgement
        // and clear latches. Only in-memory receipt checks may follow it; task
        // inputs are authenticated above and frozen by the launch authority.
        if !output.status.success() {
            return Err(io::Error::other(format!(
                "guard admission helper failed ({}); see inherited helper stderr",
                output.status
            )));
        }
        if output.stdout.len() > MAX_METADATA_BYTES as usize {
            return Err(io::Error::other(
                "guard admission receipt exceeds metadata bound",
            ));
        }
        // from_slice refuses a second JSON value or other trailing output.
        let receipt = serde_json::from_slice(&output.stdout).map_err(io::Error::other)?;
        self.validate_receipt(&receipt, worker)?;
        Ok(receipt)
    }

    fn validate_receipt(&self, receipt: &Value, worker: &WorkerIdentity) -> io::Result<()> {
        if receipt["schema"] != "vakint.resource-guard-admission.v1"
            || receipt["admitted"] != true
            || receipt["worker"]["pid"] != worker.pid
            || receipt["worker"]["start_ticks"] != worker.start_ticks
            || receipt["worker"]["uid"] != worker.uid
            || receipt["worker"]["executable_sha256"] != worker.executable_sha256
            || receipt["launch_spec_sha256"] != self.launch_spec_sha256
            || !receipt["guard"]["pid"].as_u64().is_some_and(|pid| pid > 0)
            || !receipt["guard"]["start_ticks"]
                .as_u64()
                .is_some_and(|ticks| ticks > 0)
            || receipt["guard"]["uid"] != worker.uid
            || !receipt["guard_config_sha256"]
                .as_str()
                .is_some_and(valid_sha256)
            || !receipt["slot"]
                .as_str()
                .is_some_and(|slot| !slot.is_empty() && slot.len() <= 256)
            || !receipt["admission_utc"]
                .as_str()
                .is_some_and(|utc| !utc.is_empty() && utc.len() <= 128)
            || !receipt["MemAvailable_bytes"].is_u64()
        {
            return Err(io::Error::other(
                "guard admission did not authenticate this exact Rust worker and launch specification",
            ));
        }
        Ok(())
    }
}

struct WorkerIdentity {
    pid: u32,
    start_ticks: u64,
    uid: u32,
    executable_sha256: String,
}

impl WorkerIdentity {
    fn current() -> io::Result<Self> {
        #[cfg(target_os = "linux")]
        {
            use std::os::unix::fs::MetadataExt;
            let stat = std::fs::read_to_string("/proc/self/stat")?;
            let start_ticks = stat
                .rsplit_once(')')
                .and_then(|(_, fields)| fields.split_whitespace().nth(19))
                .and_then(|ticks| ticks.parse().ok())
                .ok_or_else(|| io::Error::other("cannot authenticate Rust worker start ticks"))?;
            Ok(Self {
                pid: std::process::id(),
                start_ticks,
                uid: std::fs::metadata("/proc/self")?.uid(),
                executable_sha256: fingerprint(Path::new("/proc/self/exe"), None)?,
            })
        }
        #[cfg(not(target_os = "linux"))]
        {
            Err(io::Error::new(
                io::ErrorKind::Unsupported,
                "task guard admission requires Linux process identity",
            ))
        }
    }
}

fn valid_sha256(value: &str) -> bool {
    value.len() == 64
        && value
            .bytes()
            .all(|byte| byte.is_ascii_digit() || (b'a'..=b'f').contains(&byte))
}

fn fingerprint(path: &Path, limit: Option<u64>) -> io::Result<String> {
    let mut options = OpenOptions::new();
    options.read(true);
    #[cfg(unix)]
    {
        use std::os::unix::fs::OpenOptionsExt;
        // /proc/self/exe is the kernel-owned process symlink; bounded task inputs
        // must instead be actual regular files, never a FIFO or substituted link.
        let mut flags = rustix::fs::OFlags::NONBLOCK;
        if limit.is_some() {
            flags |= rustix::fs::OFlags::NOFOLLOW;
        }
        options.custom_flags(flags.bits() as i32);
    }
    let mut file = options.open(path)?;
    let metadata = file.metadata()?;
    if !metadata.is_file() || limit.is_some_and(|limit| metadata.len() > limit) {
        return Err(io::Error::other(
            "guard input must be a bounded regular file",
        ));
    }
    let mut digest = Sha256::new();
    let mut total = 0_u64;
    let mut buffer = [0_u8; 8192];
    loop {
        let read = file.read(&mut buffer)?;
        if read == 0 {
            break;
        }
        total = total
            .checked_add(read as u64)
            .ok_or_else(|| io::Error::other("guard file size overflow"))?;
        if limit.is_some_and(|limit| total > limit) {
            return Err(io::Error::other("guard input grew beyond metadata bound"));
        }
        digest.update(&buffer[..read]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn request() -> Request {
        Request {
            python: "/workspace/python".into(),
            helper: "/workspace/helper.py".into(),
            helper_sha256: "a".repeat(64),
            launch_spec: "/workspace/launch.json".into(),
            launch_spec_sha256: "b".repeat(64),
        }
    }

    #[test]
    fn absent_tuple_is_inert_and_every_partial_tuple_is_refused() {
        assert!(Request::from_values(Default::default()).unwrap().is_none());
        for index in 0..5 {
            let mut values: [Option<OsString>; 5] = Default::default();
            values[index] = Some("present".into());
            assert!(Request::from_values(values).is_err());
        }
        assert!(!valid_sha256("missing"));
        assert!(!valid_sha256(&"A".repeat(64)));
        let values = [
            "/workspace/python".into(),
            "/workspace/helper.py".into(),
            "a".repeat(64),
            "/workspace/launch.json".into(),
            "b".repeat(64),
        ];
        assert!(
            Request::from_values(values.clone().map(|value| Some(value.into())))
                .unwrap()
                .is_some()
        );
        let mut relative = values.clone();
        relative[1] = "relative-helper.py".into();
        assert!(Request::from_values(relative.map(|value| Some(value.into()))).is_err());
        let mut invalid = values;
        invalid[4] = "not-a-hash".into();
        assert!(Request::from_values(invalid.map(|value| Some(value.into()))).is_err());
    }

    #[test]
    fn receipt_must_match_rust_pid_start_uid_executable_and_launch_spec() {
        let request = request();
        let worker = WorkerIdentity {
            pid: 31,
            start_ticks: 42,
            uid: 7,
            executable_sha256: "c".repeat(64),
        };
        let receipt = serde_json::json!({"schema":"vakint.resource-guard-admission.v1","admitted":true,
            "worker":{"pid":31,"start_ticks":42,"uid":7,"executable_sha256":"c".repeat(64)},
            "guard":{"pid":99,"start_ticks":100,"uid":7},"launch_spec_sha256":"b".repeat(64),
            "guard_config_sha256":"d".repeat(64),"slot":"no-native-test","MemAvailable_bytes":60000000000_u64,"admission_utc":"fixture"});
        request.validate_receipt(&receipt, &worker).unwrap();
        for field in ["pid", "start_ticks", "uid", "executable_sha256"] {
            let mut bad = receipt.clone();
            bad["worker"][field] = Value::Null;
            assert!(request.validate_receipt(&bad, &worker).is_err());
        }
        for field in [
            "schema",
            "admitted",
            "launch_spec_sha256",
            "guard",
            "guard_config_sha256",
            "slot",
            "MemAvailable_bytes",
            "admission_utc",
        ] {
            let mut bad = receipt.clone();
            bad.as_object_mut().unwrap().remove(field);
            assert!(request.validate_receipt(&bad, &worker).is_err());
        }
        assert!(serde_json::from_slice::<Value>(b"{} {} ").is_err());
    }

    #[test]
    fn changed_inputs_and_nonzero_helper_never_produce_admission() {
        let directory = tempfile::tempdir().unwrap();
        let helper = directory.path().join("helper.py");
        let spec = directory.path().join("launch.json");
        std::fs::write(&helper, b"not executed by a shell").unwrap();
        std::fs::write(&spec, b"{}").unwrap();
        let mut request = Request {
            python: std::env::current_exe().unwrap(),
            helper: helper.clone(),
            helper_sha256: fingerprint(&helper, Some(MAX_METADATA_BYTES)).unwrap(),
            launch_spec: spec.clone(),
            launch_spec_sha256: fingerprint(&spec, Some(MAX_METADATA_BYTES)).unwrap(),
        };
        request.verify_inputs().unwrap();
        std::fs::write(&spec, b"changed").unwrap();
        assert!(request.verify_inputs().is_err());
        std::fs::write(&spec, b"{}").unwrap();
        request.helper_sha256 = "0".repeat(64);
        assert!(request.verify_inputs().is_err());
        request.helper_sha256 = fingerprint(&helper, Some(MAX_METADATA_BYTES)).unwrap();
        // The Rust test harness rejects Python's -B flag before running tests.
        // This is a real nonzero child, not Python/Native/guard execution.
        let worker = WorkerIdentity {
            pid: 31,
            start_ticks: 42,
            uid: 7,
            executable_sha256: "c".repeat(64),
        };
        assert!(
            request
                .admit(&worker)
                .unwrap_err()
                .to_string()
                .contains("helper failed")
        );
    }

    #[cfg(unix)]
    #[test]
    fn metadata_fifo_and_oversize_are_refused_before_hashing() {
        let directory = tempfile::tempdir().unwrap();
        let fifo = directory.path().join("helper.fifo");
        rustix::fs::mkfifoat(
            rustix::fs::CWD,
            &fifo,
            rustix::fs::Mode::RUSR | rustix::fs::Mode::WUSR,
        )
        .unwrap();
        let started = std::time::Instant::now();
        assert!(fingerprint(&fifo, Some(MAX_METADATA_BYTES)).is_err());
        assert!(started.elapsed() < std::time::Duration::from_secs(1));
        let large = directory.path().join("large.json");
        std::fs::write(&large, vec![b' '; MAX_METADATA_BYTES as usize + 1]).unwrap();
        assert!(fingerprint(&large, Some(MAX_METADATA_BYTES)).is_err());
    }

    #[test]
    #[ignore = "requires reviewed task ENV, authoritative postspawn ticket and real guard; no Native work"]
    fn real_guard_admission_returns_before_any_native_work() {
        let receipt = admit_from_task_env()
            .unwrap()
            .expect("explicit reviewed admission tuple required");
        assert_eq!(receipt["admitted"], true);
        println!(
            "{}",
            serde_json::json!({"event":"guard_admission_without_native","receipt":receipt})
        );
    }
}
