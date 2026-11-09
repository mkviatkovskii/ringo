use std::io::Write;
use std::process::{Command, Stdio};

#[allow(clippy::unwrap_used, reason = "test helper failures should panic")]
fn run(input: &str, args: &[&str]) -> std::process::Output {
    let mut child = Command::new(env!("CARGO_BIN_EXE_ringo-canonical"))
        .arg("-")
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child
        .stdin
        .take()
        .unwrap()
        .write_all(input.as_bytes())
        .unwrap();
    child.wait_with_output().unwrap()
}

#[test]
fn canonicalizes_records_from_stdin() {
    let output = run("CCO first\n\nOCC second\nC(O)C third\n", &[]);
    assert!(output.status.success(), "{output:?}");
    let text = String::from_utf8(output.stdout).unwrap();
    let lines: Vec<_> = text.lines().collect();
    assert_eq!(lines.len(), 3);
    assert!(lines.iter().all(|line| line == &lines[0]));
}

#[test]
fn invalid_record_reports_line_and_fails() {
    let output = run("CCO\n\nC(\n", &[]);
    assert!(!output.status.success());
    assert!(String::from_utf8(output.stderr).unwrap().contains("line 3"));
}

#[test]
fn exhausted_budget_fails_without_a_provisional_result() {
    let output = run("C1CCCCC1\n", &["--search-limit", "1"]);
    assert!(!output.status.success());
    assert_eq!(output.stdout, Vec::<u8>::new());
    assert!(String::from_utf8(output.stderr)
        .unwrap()
        .contains("SearchLimit"));
}
