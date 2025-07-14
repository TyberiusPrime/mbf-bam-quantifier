use std::fs::{self, File};
use std::io::Write;
use std::path::Path;
use walkdir::WalkDir;

fn main() {
    let mut out = String::new();
    out += r#"
    mod test_runner;
    use test_runner::run_test;
    "#;

    for entry in WalkDir::new("test_cases")
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| e.file_name() == "input.toml")
    {
        let case_dir = entry.path().parent().unwrap();
        if case_dir.file_name().unwrap() == "actual" {
            continue; // Skip the compare to this "actual" directory
        }
        let name = case_dir
            .strip_prefix("test_cases")
            .unwrap()
            .to_string_lossy()
            .replace(|c: char| !c.is_ascii_alphanumeric(), "_")
            .to_lowercase();

        let case_path = case_dir.to_string_lossy();

        out += &format!(
            r#"
#[test]
fn test_case_{name}() {{
    run_test(std::path::Path::new("{case_path}"));
}}
"#
        );
    }

    let out_path = Path::new("tests/generated_by_build.rs");
    fs::create_dir_all(out_path.parent().unwrap()).unwrap();
    File::create(out_path)
        .unwrap()
        .write_all(out.as_bytes())
        .unwrap();
    println!("cargo:rerun-if-changed=test_cases/");
}
