pub trait LogMessage {
    fn log_display(&self) -> String;
    fn log_file(&self) -> String {
        self.log_display()
    }
}

#[cfg(feature = "symbolica")]
mod symbolica_impl {
    use super::*;
    use spenso::shadowing::symbolica_utils::SpensoPrintSettings;
    use symbolica::prelude::*;

    impl LogMessage for Atom {
        fn log_display(&self) -> String {
            let mut settings = SpensoPrintSettings::compact().nice_symbolica();
            settings.max_line_length = Some(120);
            format!(
                "n-terms: {}, byte size: {},expr:{}",
                self.as_view().nterms(),
                self.as_view().get_byte_size(),
                self.printer(settings)
            )
        }

        fn log_file(&self) -> String {
            format!(
                "n-terms: {}, byte size: {},expr:{}",
                self.as_view().nterms(),
                self.as_view().get_byte_size(),
                self.to_plain_string()
            )
        }
    }
}
