use std::{
    fs::{File, OpenOptions},
    io::Write,
    path::PathBuf,
    sync::{Arc, Mutex, OnceLock},
    time::SystemTime,
};

use serde::Serialize;

pub static DEBUG_LOGGER: DebugLogger = DebugLogger::init();

pub struct DebugLogger {
    logger: OnceLock<File>,
}

impl DebugLogger {
    const fn init() -> Self {
        Self {
            logger: OnceLock::new(),
        }
    }

    pub fn set(&self, path: &PathBuf) -> Result<(), File> {
        let file = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)
            .unwrap();

        self.logger.set(file)
    }

    #[cold]
    pub fn write<T: Serialize>(&self, msg: &str, data: &T) {
        let log_message = LogMessage {
            msg,
            data,
            ts: SystemTime::now(),
        };

        let json_log_message = serde_json::to_string(&log_message).unwrap();
        writeln!(self.logger.get().unwrap(), "{}", json_log_message);
    }
}

#[derive(Serialize, Clone, Debug)]
struct LogMessage<'a, T> {
    msg: &'a str,
    ts: std::time::SystemTime,
    data: &'a T,
}
