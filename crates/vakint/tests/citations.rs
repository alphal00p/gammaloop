#![cfg(all(unix, feature = "symbolica_community_module"))]

use symbolica::api::python::SymbolicaCommunityModule;
use vakint::{Vakint, VakintSettings, symbolica_community_module::VakintWrapper};

#[test]
fn backend_execution_updates_community_citations() {
    assert!(VakintWrapper::get_citations().is_empty());
    let vakint = Vakint::new().unwrap();
    let settings = VakintSettings {
        python_exe_path: "sh".into(),
        form_exe_path: "sh".into(),
        ..VakintSettings::default()
    };
    assert_eq!(VakintWrapper::get_citations().len(), 1);

    // Exercise the real subprocess adapters without numerical pySecDec or FORM.
    // A second call uses another Vakint instance to check process-wide tracking.
    let input = (
        "run.sh".into(),
        "printf 'adapter result' > out.txt\n".into(),
    );
    for engine in [&vakint, &vakint.clone()] {
        let result = engine
            .run_pysecdec(
                &settings,
                std::slice::from_ref(&input),
                vec![],
                true,
                None,
                None,
            )
            .unwrap();
        assert_eq!(result, ["adapter result"]);
        let citations = VakintWrapper::get_citations();
        assert_eq!(citations.len(), 2);
        assert_eq!(citations[1].id, "arXiv:1703.09692");
    }

    let result = vakint
        .run_form(&settings, &[], input, vec![], true, None)
        .unwrap();
    assert_eq!(result, "adapter result");
    let citations = VakintWrapper::get_citations();
    assert_eq!(citations.len(), 3);
    assert_eq!(citations[1].id, "arXiv:1203.6543");
    assert_eq!(citations[2].id, "arXiv:1703.09692");
}
