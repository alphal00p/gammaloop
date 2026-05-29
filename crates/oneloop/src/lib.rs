//! `oneloop`: one-loop IBP reduction

    env!("CARGO_PKG_VERSION")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn builds_and_reports_version() {
        assert_eq!(version(), "0.1.0");
    }
}
