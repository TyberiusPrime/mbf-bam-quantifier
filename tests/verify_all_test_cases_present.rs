use std::collections::HashSet;
use std::fs;
use walkdir::WalkDir;

#[test]
fn all_test_cases_are_generated() {
    let generated = fs::read_to_string("tests/generated.rs")
        .expect("Failed to read generated.rs");

    let mut expected_tests = HashSet::new();

    for entry in WalkDir::new("test_cases")
        .into_iter()
        .filter_map(Result::ok)
        .filter(|e| e.file_name() == "input.toml")
    {
        let case_dir = entry.path().parent().unwrap();
        if case_dir.file_name().unwrap() == "actual" {
            continue;
        }

        let name = case_dir
            .strip_prefix("test_cases")
            .unwrap()
            .to_string_lossy()
            .replace(|c: char| !c.is_ascii_alphanumeric(), "_")
            .to_lowercase();

        expected_tests.insert(format!("fn test_case_{}()", name));
    }

    for test_fn in expected_tests {
        assert!(
            generated.contains(&test_fn),
            "Missing test function: {}. Rerun ./dev/update_tests.py",
            test_fn
        );
    }
}

