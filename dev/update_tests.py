#!/usr/bin/env python3

from pathlib import Path


assert Path("test_cases").exists(), "Starting from the wrong dir, test_cases not found"

out = """
mod test_runner;
use test_runner::run_test;
"""
for input_toml in sorted(Path("test_cases").rglob("input.toml")):
    case_dir = Path(input_toml.parent)
    if case_dir.name == "actual":
        continue

    name = str(case_dir.relative_to("test_cases"))
    name = "".join(c if c.isascii() and c.isalnum() else "_" for c in name).lower()
    case_path = str(case_dir)

    out += f'''
#[test]
fn test_case_{name}() {{
    run_test(std::path::Path::new("{case_path}"));
}}
'''

out_path = Path("tests/generated.rs")
out_path.parent.mkdir(parents=True, exist_ok=True)
out_path.write_text(out)
