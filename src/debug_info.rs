use std::{
    fs::OpenOptions,
    path::PathBuf,
    sync::{Arc, Mutex, OnceLock},
};

use slog::{o, Drain, Logger};

pub static DEBUG_LOGGER: DebugLogger = DebugLogger::init();

pub struct DebugLogger {
    logger: OnceLock<Logger>,
}

impl DebugLogger {
    const fn init() -> Self {
        Self {
            logger: OnceLock::new(),
        }
    }

    pub fn set(&self, path: &PathBuf) -> Result<(), Logger> {
        let file = OpenOptions::new()
            .create(true)
            .write(true)
            .truncate(true)
            .open(path)
            .unwrap();

        let drain = Arc::new(Mutex::new(slog_json::Json::default(file))).fuse();
        let logger = slog::Logger::root(drain, o!());
        self.logger.set(logger)
    }

    #[cold]
    pub fn get(&self) -> &Logger {
        self.logger.get().unwrap()
    }
}
