use std::{
    collections::hash_map::Entry::Vacant,
    ffi::OsStr,
    fs::{create_dir_all, File, OpenOptions},
    io::Write,
    path::PathBuf,
    sync::{LazyLock, Mutex},
    time::SystemTime,
};

use ahash::HashMap;
use color_eyre::Report;
use eyre::eyre;
use serde::Serialize;

use crate::{Precision, RotationMethod};

pub static DEBUG_LOGGER: DebugLogger = DebugLogger::init();

/// This could also include channel_id or graph_id if the user is doing the explicit sum over lmb channels/graphs, but I don't
/// support it for now, as it is not used much especially not when debugging.
#[derive(Hash, Copy, Clone, Debug, Eq, PartialEq)]
pub enum EvalState {
    General,
    PrecRot((RotationMethod, Precision)),
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

impl DebugLogger {
    const fn init() -> Self {
        Self {
            logger: LazyLock::new(|| {
                Mutex::new(LogImpl {
                    files: HashMap::default(),
                    current: EvalState::General,
                    path: PathBuf::new(),
                })
            }),
        }
    }

    #[cold]
    pub fn new_log(&self, path: &PathBuf) -> Result<(), Report> {
        self.logger.lock().unwrap().new_log(path)
    }

    #[cold]
    pub fn write<T: Serialize>(&self, msg: &str, data: &T) {
        self.logger
            .lock()
            .unwrap()
            .write(msg, data)
            .expect("failed to write log")
    }

    #[cold]
    pub fn set_state(&self, state: EvalState) {
        self.logger
            .lock()
            .unwrap()
            .set_state(state)
            .expect("failed to change state")
    }
}

struct LogImpl {
    files: HashMap<EvalState, File>,
    current: EvalState,
    path: PathBuf,
}

impl LogImpl {
    fn new_log(&mut self, path: &PathBuf) -> Result<(), Report> {
        self.files.clear();
        self.current = EvalState::General;
        self.path = path.clone();

        let extension = path.extension().and_then(OsStr::to_str).map_or_else(
            || Err(eyre!("could not determine extension, is it a .glog?")),
            Ok,
        )?;

        assert_eq!(extension, "glog");

        create_dir_all(path)?;
        let general_path = path.join("general.jsonl");

        let file = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(general_path)?;

        self.files.insert(EvalState::General, file);

        Ok(())
    }

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
}

#[derive(Serialize, Clone, Debug)]
struct LogMessage<'a, T> {
    msg: &'a str,
    ts: std::time::SystemTime,
    data: &'a T,
}
