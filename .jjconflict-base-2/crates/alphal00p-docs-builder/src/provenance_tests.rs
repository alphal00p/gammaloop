use std::{env, fs, process::Command};

use super::SiteBuilder;

#[test]
fn provenance_uses_the_documented_checkout_and_explicit_overrides() {
    const CASE: &str = "ALPHAL00P_DOCS_PROVENANCE_TEST_CASE";
    let Ok(case) = env::var(CASE) else {
        // Isolate process-wide overrides without mutating the parallel test runner's environment.
        for case in ["checkout", "explicit", "ci"] {
            let mut command = Command::new(env::current_exe().unwrap());
            command
                .args([
                    "--exact",
                    "provenance_tests::provenance_uses_the_documented_checkout_and_explicit_overrides",
                    "--nocapture",
                ])
                .env(CASE, case);
            for variable in [
                "ALPHAL00P_DOCS_GIT_COMMIT",
                "GITHUB_SHA",
                "ALPHAL00P_DOCS_GIT_TIMESTAMP",
                "SOURCE_DATE_EPOCH",
                "GIT_DIR",
                "GIT_WORK_TREE",
            ] {
                command.env_remove(variable);
            }
            if case != "checkout" {
                command
                    .env("GITHUB_SHA", "ci-commit")
                    .env("SOURCE_DATE_EPOCH", "123");
            }
            if case == "explicit" {
                command
                    .env("ALPHAL00P_DOCS_GIT_COMMIT", " explicit-commit ")
                    .env("ALPHAL00P_DOCS_GIT_TIMESTAMP", " 456 ");
            }
            assert!(command.status().unwrap().success(), "case: {case}");
        }
        return;
    };

    let temporary = tempfile::tempdir().unwrap();
    let mut builder = SiteBuilder::discover().unwrap();
    builder.root = temporary.path().to_path_buf();
    if case != "checkout" {
        // Exported sources have no repository metadata; explicit values remain sufficient.
        let expected = if case == "explicit" {
            ("explicit-commit", 456)
        } else {
            ("ci-commit", 123)
        };
        assert_eq!(builder.git_commit().unwrap(), expected.0);
        assert_eq!(builder.git_timestamp().unwrap(), expected.1);
        return;
    }

    let error = builder.git_commit().unwrap_err().to_string();
    assert!(error.contains("git rev-parse HEAD"), "{error}");
    assert!(error.contains("ALPHAL00P_DOCS_GIT_COMMIT"), "{error}");
    let error = builder.git_timestamp().unwrap_err().to_string();
    assert!(error.contains("SOURCE_DATE_EPOCH"), "{error}");

    for kind in ["git", "standalone", "linked", "colocated"] {
        builder.root = temporary.path().join(kind);
        fs::create_dir(&builder.root).unwrap();
        let initialization: &[&[&str]] = match kind {
            "git" => &[
                &["git", "init"],
                &[
                    "git",
                    "-c",
                    "user.name=Docs Test",
                    "-c",
                    "user.email=docs@example.invalid",
                    "commit",
                    "--allow-empty",
                    "-m",
                    "initial",
                ],
            ],
            "standalone" => &[&["jj", "git", "init", "--no-colocate"]],
            "linked" => &[&[
                "jj",
                "-R",
                "../standalone",
                "workspace",
                "add",
                "--name",
                "linked",
                ".",
            ]],
            "colocated" => &[&["jj", "git", "init", "--colocate"]],
            _ => unreachable!(),
        };
        for args in initialization {
            let output = Command::new(args[0])
                .args(&args[1..])
                .current_dir(&builder.root)
                .env("GIT_AUTHOR_DATE", "2025-01-02T03:04:05+00:00")
                .env("GIT_COMMITTER_DATE", "2025-01-02T03:04:05+00:00")
                .output()
                .unwrap();
            assert!(output.status.success(), "{kind}: {args:?}: {output:?}");
        }
        if kind != "git" {
            let output = Command::new("jj")
                .args([
                    "--config",
                    "user.name=Docs Test",
                    "--config",
                    "user.email=docs@example.invalid",
                    "--config",
                    "debug.commit-timestamp=2025-01-02T03:04:05+00:00",
                    "describe",
                    "-m",
                    "documented change",
                ])
                .current_dir(&builder.root)
                .output()
                .unwrap();
            assert!(output.status.success(), "{kind}: {output:?}");
        }
        if kind == "colocated" {
            // Give Git HEAD a real parent revision that differs from the documented JJ change.
            let output = Command::new("jj")
                .args([
                    "--config",
                    "user.name=Docs Test",
                    "--config",
                    "user.email=docs@example.invalid",
                    "--config",
                    "debug.commit-timestamp=2025-01-02T03:04:05+00:00",
                    "new",
                    "-m",
                    "current documentation",
                ])
                .current_dir(&builder.root)
                .output()
                .unwrap();
            assert!(output.status.success(), "{output:?}");
        }
        if matches!(kind, "standalone" | "linked") {
            assert!(!builder.root.join(".git").exists(), "{kind}");
        }
        if kind == "linked" {
            assert!(builder.root.join(".jj/repo").is_file());
        }
        let args: &[&str] = if kind == "git" {
            &["git", "rev-parse", "HEAD"]
        } else {
            &[
                "jj",
                "--ignore-working-copy",
                "log",
                "-r",
                "@",
                "--no-graph",
                "-T",
                "commit_id",
            ]
        };
        let output = Command::new(args[0])
            .args(&args[1..])
            .current_dir(&builder.root)
            .output()
            .unwrap();
        assert!(output.status.success(), "{kind}: {output:?}");
        let expected_commit = String::from_utf8(output.stdout).unwrap().trim().to_owned();
        assert_eq!(expected_commit.len(), 40, "{kind}");
        fs::write(builder.root.join("uncommitted.txt"), "work in progress").unwrap();
        // Reading provenance must not snapshot edits or rewrite the working-copy commit.
        assert_eq!(builder.git_commit().unwrap(), expected_commit, "{kind}");
        assert_eq!(builder.git_timestamp().unwrap(), 1_735_787_045, "{kind}");
        if kind == "colocated" {
            let output = Command::new("git")
                .args(["rev-parse", "HEAD"])
                .current_dir(&builder.root)
                .output()
                .unwrap();
            assert!(output.status.success(), "{output:?}");
            assert_ne!(
                String::from_utf8(output.stdout).unwrap().trim(),
                expected_commit
            );
        }
    }

    builder.root = temporary.path().join("broken-jj");
    fs::create_dir_all(builder.root.join(".jj")).unwrap();
    let error = builder.git_commit().unwrap_err().to_string();
    assert!(error.contains("jj"), "{error}");
    assert!(
        error.contains(&builder.root.display().to_string()),
        "{error}"
    );
}

#[test]
#[cfg(unix)]
fn each_build_captures_one_revision_before_rendering() {
    use std::os::unix::{fs::PermissionsExt, fs::symlink};

    use super::{BuildChannel, BuildRequest, RustdocCacheMode};

    const CASE: &str = "ALPHAL00P_DOCS_BUILD_PROVENANCE_TEST_ROOT";
    let Ok(root) = env::var(CASE) else {
        let temporary = tempfile::tempdir().unwrap();
        let checkout = SiteBuilder::discover().unwrap();
        let root = temporary.path().join("checkout");
        fs::create_dir(&root).unwrap();
        for entry in fs::read_dir(checkout.root).unwrap() {
            let entry = entry.unwrap();
            if !matches!(entry.file_name().to_str(), Some(".git" | ".jj" | "target")) {
                symlink(entry.path(), root.join(entry.file_name())).unwrap();
            }
        }
        fs::create_dir(root.join(".jj")).unwrap();
        let bin = temporary.path().join("bin");
        fs::create_dir(&bin).unwrap();
        let shell = Command::new("sh")
            .args(["-c", "command -v sh"])
            .output()
            .unwrap();
        assert!(shell.status.success());
        let shell = String::from_utf8(shell.stdout).unwrap();
        let executable = bin.join("jj");
        fs::write(
            &executable,
            format!(
                r#"#!{}
case "$*" in
    *committer.timestamp*)
        case "$*" in
            *1111111111111111111111111111111111111111*) printf 123 ;;
            *2222222222222222222222222222222222222222*) printf 456 ;;
            *) printf 'timestamp must use the captured revision' >&2; exit 1 ;;
        esac ;;
    *)
        if [ -f "$ALPHAL00P_DOCS_BUILD_PROVENANCE_TEST_ROOT/read-commit" ]; then
            printf 2222222222222222222222222222222222222222
        else
            : > "$ALPHAL00P_DOCS_BUILD_PROVENANCE_TEST_ROOT/read-commit"
            printf 1111111111111111111111111111111111111111
        fi ;;
esac
"#,
                shell.trim()
            ),
        )
        .unwrap();
        fs::set_permissions(executable, fs::Permissions::from_mode(0o755)).unwrap();
        let inherited_path = env::var_os("PATH").unwrap();
        let path =
            env::join_paths(std::iter::once(bin).chain(env::split_paths(&inherited_path))).unwrap();
        let mut command = Command::new(env::current_exe().unwrap());
        command
            .args([
                "--exact",
                "provenance_tests::each_build_captures_one_revision_before_rendering",
                "--nocapture",
            ])
            .env(CASE, &root)
            .env("PATH", path);
        for variable in [
            "ALPHAL00P_DOCS_GIT_COMMIT",
            "GITHUB_SHA",
            "ALPHAL00P_DOCS_GIT_TIMESTAMP",
            "SOURCE_DATE_EPOCH",
        ] {
            command.env_remove(variable);
        }
        assert!(command.status().unwrap().success());
        return;
    };

    let mut builder = SiteBuilder::discover().unwrap();
    builder.root = root.into();
    for (commit, timestamp) in [
        ("1111111111111111111111111111111111111111", 123),
        ("2222222222222222222222222222222222222222", 456),
    ] {
        let output = builder.root.join("target/site");
        builder
            .build(BuildRequest {
                product: "linnet".to_owned(),
                channel: BuildChannel::Latest,
                output: output.clone(),
                snapshot_tag: None,
                include_rustdoc: false,
                include_typst: true,
                rustdoc_target_root: None,
                rustdoc_cache: RustdocCacheMode::Disabled,
                dependency_output: None,
            })
            .unwrap();
        let metadata: serde_json::Value = serde_json::from_slice(
            &fs::read(output.join("products/linnet/latest/snapshot.json")).unwrap(),
        )
        .unwrap();
        assert_eq!(metadata["git_commit"], commit);
        assert_eq!(metadata["git_timestamp"], timestamp);
    }
}
