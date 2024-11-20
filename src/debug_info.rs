use std::{
    collections::hash_map::Entry::Vacant,
    fmt::Debug,
    fs::{self, create_dir_all, File, OpenOptions},
    hash::Hash,
    io::Write,
    path::PathBuf,
    sync::{LazyLock, Mutex},
    time::SystemTime,
};

use ahash::HashMap;
use color_eyre::Report;
use serde::Serialize;

use crate::{Precision, RotationSetting};

pub static DEBUG_LOGGER: DebugLogger = DebugLogger::init();

/// This could also include channel_id or graph_id if the user is doing the explicit sum over lmb channels/graphs, but I don't
/// support it for now, as it is not used much especially not when debugging.
#[derive(Copy, Clone, Debug)]
pub enum EvalState {
    General,
    PrecRot((RotationSetting, Precision)),
}

impl From<&EvalState> for String {
    fn from(value: &EvalState) -> Self {
        match value {
            EvalState::General => "general".to_owned(),
            EvalState::PrecRot((rotation_settings, prec)) => {
                format!("{}_{}", rotation_settings.as_str(), prec.as_string())
            }
        }
    }
}

// This is needed since RotationSetting may now contain f64
impl PartialEq for EvalState {
    fn eq(&self, other: &Self) -> bool {
        Into::<String>::into(self) == Into::<String>::into(other)
    }
}

impl Eq for EvalState {}

impl Hash for EvalState {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        String::from(self).hash(state);
    }
}

impl EvalState {
    fn as_file_path(&self) -> PathBuf {
        match self {
            Self::General => PathBuf::from("general.jsonl"),
            Self::PrecRot((rotation_method, prec)) => {
                PathBuf::from(format!("{}_{}.jsonl", rotation_method.as_str(), prec,))
            }
        }
    }
}

pub struct DebugLogger {
    logger: LazyLock<Mutex<LogImpl>>,
}

impl momtrop::log::Logger for DebugLogger {
    fn write<T: Serialize>(&self, msg: &str, data: &T) {
        self.write(msg, data)
    }
}

impl DebugLogger {
    const fn init() -> Self {
        Self {
            logger: LazyLock::new(|| {
                let log_impl = LogImpl::new();
                Mutex::new(log_impl)
            }),
        }
    }

    #[cold]
    pub fn write<T: Serialize>(&self, msg: &str, data: &T) {
        self.logger
            .lock()
            .unwrap()
            .write(msg, data)
            .unwrap_or_else(|_| panic!("failed to write log for log message: {}", msg))
    }

    #[cold]
    pub fn set_state(&self, state: EvalState) {
        self.logger
            .lock()
            .unwrap()
            .set_state(state)
            .expect("failed to change state")
    }

    pub fn reset(&self) {
        self.logger.lock().unwrap().reset();
    }
}

struct LogImpl {
    files: HashMap<EvalState, File>,
    current: EvalState,
    path: PathBuf,
}

impl LogImpl {
    fn set_state(&mut self, state: EvalState) -> Result<(), Report> {
        self.current = state;
        if let Vacant(e) = self.files.entry(state) {
            let path = self.path.join(state.as_file_path());

            let file = OpenOptions::new()
                .create(true)
                .write(true)
                .truncate(true)
                .open(path)?;

            e.insert(file);

            Ok(())
        } else {
            Ok(())
        }
    }

    fn write<T: Serialize>(&mut self, msg: &str, data: &T) -> Result<(), Report> {
        let log_message = LogMessage {
            msg,
            data,
            ts: SystemTime::now(),
        };

        let json_log_message = serde_json::to_string(&log_message)?;

        writeln!(
            self.files.get(&self.current).unwrap(),
            "{}",
            json_log_message
        )?;

        Ok(())
    }
    fn reset(&mut self) {
        *self = Self::new();
    }

    fn new() -> Self {
        let path = PathBuf::from("log.glog");
        let path_to_general = path.join("general.jsonl");

        if path.exists() {
            fs::remove_dir_all(&path).expect("unable to remove older log file");
        }

        create_dir_all(&path).expect("unable to create log file");

        let general_file = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path_to_general)
            .expect("unable to create log file");

        let mut files = HashMap::default();
        files.insert(EvalState::General, general_file);

        LogImpl {
            files,
            current: EvalState::General,
            path,
        }
    }
}

#[derive(Serialize, Clone, Debug)]
struct LogMessage<'a, T> {
    msg: &'a str,
    ts: std::time::SystemTime,
    data: &'a T,
}
