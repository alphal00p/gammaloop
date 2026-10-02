//! Exercise registration before any tensor or symbolic state is warmed.

use std::{
    process::Command,
    thread,
    time::{Duration, Instant},
};

#[test]
fn vector_registration_bootstraps_in_a_fresh_process() {
    const CHILD: &str = "FEYNKIT_VECTOR_BOOTSTRAP_CHILD";
    if std::env::var_os(CHILD).is_some() {
        // This must be the first Symbolica operation in this process. Linking
        // the tensor reducer also registers its dependency initializers.
        let vector = spenso::vector_symbol!("vector_bootstrap::first");
        let repeated = spenso::vector_symbol!("vector_bootstrap::first");
        assert_eq!(vector, repeated);
        assert!(vector.has_tag(&spenso::network::tags::SPENSO_TAG.rank1));
        let _ = feynkit_tensor::TensorReducer::new(symbolica::atom::Atom::num(4));
        return;
    }

    let mut child = Command::new(std::env::current_exe().unwrap())
        .args([
            "--exact",
            "vector_registration_bootstraps_in_a_fresh_process",
        ])
        .env(CHILD, "1")
        .spawn()
        .unwrap();
    let deadline = Instant::now() + Duration::from_secs(15);
    loop {
        if let Some(status) = child.try_wait().unwrap() {
            assert!(
                status.success(),
                "fresh registration process failed: {status}"
            );
            break;
        }
        if Instant::now() >= deadline {
            child.kill().unwrap();
            child.wait().unwrap();
            panic!("first vector registration deadlocked during Symbolica bootstrap");
        }
        thread::sleep(Duration::from_millis(20));
    }
}
