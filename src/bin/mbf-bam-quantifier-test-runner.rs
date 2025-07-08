use anyhow::{bail, Context, Result};
use std::fs::{self, DirEntry};
use std::io::Read;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process;
use tempfile::TempDir;

const CLI_UNDER_TEST: &str = "mbf-bam-quantifier";

fn main() -> Result<()> {
    //human_panic::setup_panic!();
    for test_dir in std::env::args().skip(1).filter(|x| !x.starts_with("--")) {
        run_tests(PathBuf::from(test_dir), false)?
    }
    if std::env::args().count() < 2 {
        let test_dir = std::env::args().nth(1).unwrap_or("test_cases".to_string());
        run_tests(PathBuf::from(test_dir), false)?
    }
    Ok(())
}

fn run_tests(test_dir: impl AsRef<Path>, continue_upon_failure: bool) -> Result<()> {
    let last_failed_filename: PathBuf =
        format!("/tmp/.{CLI_UNDER_TEST}-test-runner-last-failed").into();
    let last_failed: Option<PathBuf> = if last_failed_filename.exists() {
        Some(
            fs::read_to_string(&last_failed_filename)
                .context("Read last failed test case")?
                .trim()
                .into(),
        )
    } else {
        None
    };
    // Find test cases
    let test_dir = test_dir.as_ref();
    let mut test_cases = discover_test_cases(test_dir)?;

    //randomize order
    use rand::seq::SliceRandom;
    let mut rng = rand::rng();
    test_cases.shuffle(&mut rng);

    if let Some(last_failed) = last_failed {
        //put last failed test to the front - if present
        if test_cases.iter().any(|x| x.dir == last_failed) {
            println!(
                "Found last failed test case: {}. Running it first.",
                last_failed.display()
            );
            test_cases.retain(|x| x.dir != last_failed);
            test_cases.insert(0, TestCase::new(last_failed));
        }
    }

    let mut passed = 0;
    let mut failed = 0;
    let processor_path = find_processor()?;
    let start = std::time::Instant::now();

    println!("Found {} test cases", test_cases.len());
    for test_case in test_cases {
        if test_case.dir.join("skip").exists() {
            println!(
                "Skipping test case: {} (skip file present)",
                test_case.dir.display()
            );
            continue;
        }

        let repeat_count = fs::read_to_string(test_case.dir.join("repeat"))
            .map(|x| {
                x.trim()
                    .parse::<usize>()
                    .expect("Repeat file with non number")
            })
            .unwrap_or(1);

        for repeat in 0..repeat_count {
            let start = std::time::Instant::now();
            let test_result = if test_case.is_panic {
                print!("\n  Running panic test: {} {}", test_case.dir.display(), repeat);
                run_panic_test(&test_case, processor_path.as_ref())
            } else {
                print!("\n  Running regular test: {} {}", test_case.dir.display(), repeat);
                run_output_test(&test_case, processor_path.as_ref())
            };
            let elapsed = start.elapsed();
            print!(" ({}.{:03}s)", elapsed.as_secs(), elapsed.subsec_millis());

            match test_result {
                Ok(()) => {
                    //put checkmark before last line written
                    //so we need minimal lines, but report what we're running
                    print!("\r✅");

                    //println!("✅ Output test passed");
                    passed += 1;
                }
                Err(e) => {
                    //write last failed to file
                    std::fs::write(
                        &last_failed_filename,
                        test_case.dir.to_string_lossy().to_string(),
                    )
                    .ok();
                    print!("\r❌");
                    print!("\n{:?}", e);
                    failed += 1;
                    break; // no more repeats for this one
                }
            }
        }
        if failed > 0 && !continue_upon_failure {
            println!(
                "Stopping due to failure in test: {}",
                test_case.dir.display()
            );
            break;
        }
    }

    let elapsed = start.elapsed();
    println!(
        "\nTest results: {} passed, {} failed. Took {}.{:03}s.",
        passed,
        failed,
        elapsed.as_secs(),
        elapsed.subsec_millis()
    );

    if failed > 0 {
        process::exit(1);
    }

    Ok(())
}
///
/// Finds the full path of a binary in $PATH
fn find_in_path(bin: &str) -> Option<PathBuf> {
    std::env::var_os("PATH")?
        .to_string_lossy()
        .split(':')
        .map(PathBuf::from)
        .find_map(|dir| {
            let full_path = dir.join(bin);
            if full_path.is_file()
                && fs::metadata(&full_path).ok()?.permissions().mode() & 0o111 != 0
            {
                Some(full_path)
            } else {
                None
            }
        })
}

fn find_processor() -> Result<PathBuf> {
    // prefer the one in path
    // if it exists, use that one
    if let Some(path) = find_in_path(CLI_UNDER_TEST) {
        return Ok(path);
    }
    // otherwise, check if we have a binary next to us
    let current_exe = std::env::current_exe().context("Get current executable path")?;
    let parent = current_exe
        .parent()
        .context("Get parent directory of executable")?;
    if parent.file_name().unwrap().to_string_lossy() == "debug" {
        // run a quick cargo build in debug mod
        std::process::Command::new("cargo")
            .arg("build")
            .status()
            .context("Failed to run cargo build")?
            .success()
            .then_some(())
            .ok_or_else(|| anyhow::anyhow!("Cargo build failed"))?;
    } else if parent.file_name().unwrap().to_string_lossy() == "release" {
        // run a quick cargo build in release mod
        std::process::Command::new("cargo")
            .arg("build")
            .arg("--release")
            .status()
            .context("Failed to run cargo build")?
            .success()
            .then_some(())
            .ok_or_else(|| anyhow::anyhow!("Cargo build failed"))?;
    }
    let bin_path = current_exe
        .parent()
        .context("Get parent directory of executable")?
        .join(CLI_UNDER_TEST);

    if !bin_path.exists() {
        anyhow::bail!(
            "{CLI_UNDER_TEST} binary not found at: {}",
            bin_path.display()
        );
    }

    Ok(bin_path)
}
struct TestCase {
    dir: PathBuf,
    is_panic: bool,
}

impl TestCase {
    fn new(dir: PathBuf) -> Self {
        let is_panic = dir.join("expected_panic.txt").exists();
        TestCase { dir, is_panic }
    }
}

fn discover_test_cases(dir: &Path) -> Result<Vec<TestCase>> {
    if !dir.exists() {
        anyhow::bail!("Test directory does not exist: {}", dir.display());
    }

    let mut test_cases = Vec::new();
    discover_test_cases_recursive(dir, &mut test_cases)?;
    Ok(test_cases)
}

fn discover_test_cases_recursive(dir: &Path, test_cases: &mut Vec<TestCase>) -> Result<()> {
    // Check if this directory is a test case
    if dir.join("input.toml").exists() && !dir.join("ignore").exists()
        || dir.file_name().unwrap().to_string_lossy() == "actual"
    {
        test_cases.push(TestCase::new(dir.to_path_buf()));
        return Ok(());
    }

    // Otherwise, search through subdirectories
    for entry in fs::read_dir(dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_dir() {
            discover_test_cases_recursive(&path, test_cases)?;
        }
    }

    Ok(())
}

fn read_compressed(filename: impl AsRef<Path>) -> Result<String> {
    let fh = std::fs::File::open(filename.as_ref())
        .with_context(|| format!("Could not open file {:?}", filename.as_ref()))?;
    let mut wrapped = niffler::send::get_reader(Box::new(fh))?;
    let mut out: Vec<u8> = Vec::new();
    wrapped.0.read_to_end(&mut out)?;
    Ok(std::str::from_utf8(&out)?.to_string())
}

struct TestOutput {
    stdout: String,
    stderr: String,
    return_code: i32,
    missing_files: Vec<String>,
    mismatched_files: Vec<(String, String)>,
    unexpected_files: Vec<String>,
}

fn run_panic_test(the_test: &TestCase, processor_cmd: &Path) -> Result<()> {
    let rr = perform_test(the_test, processor_cmd)?;
    if rr.return_code == 0 {
        bail!("No panic occurred, but expected one.");
    }
    let expected_panic_file = the_test.dir.join("expected_panic.txt");
    let expected_panic_content = fs::read_to_string(&expected_panic_file)
        .context("Read expected panic file")?
        .trim()
        .to_string();

    if !rr.stderr.contains(&expected_panic_content) {
        anyhow::bail!(
            "{CLI_UNDER_TEST} did not panic as expected.\nExpected panic: {}\nActual stderr: '{}'",
            expected_panic_content,
            rr.stderr
        );
    }
    Ok(())
}

fn run_output_test(test_case: &TestCase, processor_cmd: &Path) -> Result<()> {
    let rr = perform_test(test_case, processor_cmd)?;

    if rr.return_code != 0 {
        anyhow::bail!(
            "{CLI_UNDER_TEST} failed with return code: {}\nstdout: {}\nstderr: {}",
            rr.return_code,
            rr.stdout,
            rr.stderr
        );
    }

    let mut msg = String::new();
    for missing_file in &rr.missing_files {
        msg.push_str(&format!(
            "\t- Expected output file not created: {}\n",
            missing_file
        ));
    }
    for unexpected_file in &rr.unexpected_files {
        msg.push_str(&format!(
            "\t- Unexpected output file created: {}\n",
            unexpected_file
        ));
    }
    for (actual_path, _expected_path) in &rr.mismatched_files {
        msg.push_str(&format!("\t- {} (mismatched)\n", actual_path));
    }
    if !msg.is_empty() {
        anyhow::bail!("\toutput files failed verification.\n{}", msg);
    }
    Ok(())
}

fn visit_dirs(dir: &Path, cb: &mut dyn FnMut(&DirEntry) -> Result<()>) -> Result<()> {
    if dir.is_dir() {
        for entry in fs::read_dir(dir)? {
            let entry = entry?;
            let path = entry.path();
            if path.is_dir() {
                visit_dirs(&path, cb)?;
            } else {
                cb(&entry)?;
            }
        }
    }
    Ok(())
}

fn scan_dir<F: Fn(&str, &str) -> bool>(dir: &Path, callback: F) -> Result<Vec<(PathBuf, String)>> {
    let mut files = Vec::new();
    visit_dirs(dir, &mut |entry: &DirEntry| -> Result<()> {
        let path = entry.path();
        let relative_path = path
            .strip_prefix(dir)
            .context("Strip prefix from directory path")?
            .to_string_lossy()
            .to_string();

        let path = if path.is_symlink() {
            let res_path = fs::read_link(&path)
                .with_context(|| format!("Failed to read symlink: {}", path.display()))?;
            if res_path.is_relative() {
                path.parent()
                    .context("Get parent directory of symlink")?
                    .join(res_path)
                    .canonicalize()
                    .with_context(|| format!("Failed to follow symlink: {}", path.display()))?
            } else {
                res_path
            }
        } else {
            path
        };

        if path.is_file() {
            if let Some(file_name) = path.file_name() {
                let file_name_str = file_name.to_string_lossy();
                if callback(relative_path.as_ref(), &file_name_str) {
                    files.push((path, relative_path));
                }
            }
        }
        Ok(())
    })?;
    Ok(files)
}

fn perform_test(test_case: &TestCase, processor_cmd: &Path) -> Result<TestOutput> {
    let mut result = TestOutput {
        stdout: String::new(),
        stderr: String::new(),
        return_code: 0,
        missing_files: Vec::new(),
        mismatched_files: Vec::new(),
        unexpected_files: Vec::new(),
    };

    let actual_dir = test_case.dir.join("actual");
    // Create actual directory and copy files
    if actual_dir.exists() {
        fs::remove_dir_all(&actual_dir)?;
    }

    let input_files = scan_dir(test_case.dir.as_path(), |relative_path, _filename| {
        relative_path.starts_with("input")
    })?;
    let expected_files = scan_dir(test_case.dir.as_path(), |relative_path, filename| {
        !filename.starts_with("input")
            && !filename.starts_with("ignore_")
            && !relative_path.starts_with("ignore_")
    })?;

    let temp_dir = setup_test_environment(input_files).context("Setup test dir")?;

    // Run the processor
    let config_file = temp_dir.path().join("input.toml");
    //chdir to temp_dir

    fs::create_dir_all(&actual_dir)?;
    //copy all files from temp_dir to actual_dir
    visit_dirs(temp_dir.path(), &mut |entry: &DirEntry| -> Result<()> {
        let path = entry.path();
        let relative_path = path
            .strip_prefix(temp_dir.path())
            .context("Strip prefix from temp dir path")?;
        let target = actual_dir.join(relative_path);
        std::fs::create_dir_all(target.parent().unwrap())
            .context("Create parent directory for actual file")?;
        if path.is_file() {
            fs::copy(&path, &target)?;
        }
        Ok(())
    })?;

    let proc = std::process::Command::new(processor_cmd)
        .arg(&config_file)
        //.arg(temp_dir.path())
        .env("NO_FRIENDLY_PANIC", "1")
        .current_dir(temp_dir.path())
        .output()
        .context(format!("Failed to run {CLI_UNDER_TEST}"))?;

    let all_files_in_temp_dir = scan_dir(temp_dir.path(), |_, _| true)?;
    copy_files(&all_files_in_temp_dir, actual_dir.as_path())?;

    let stdout = String::from_utf8_lossy(&proc.stdout);
    let stderr = String::from_utf8_lossy(&proc.stderr);
    result.return_code = proc.status.code().unwrap_or(-1);
    result.stdout = stdout.to_string();
    result.stderr = stderr.to_string();

    //for debugging..
    fs::write(actual_dir.as_path().join("stdout"), stdout.as_bytes())
        .context("Failed to write stdout to file")?;
    fs::write(actual_dir.as_path().join("stderr"), stderr.as_bytes())
        .context("Failed to write stderr to file")?;

    let output_files_in_temp_dir = scan_dir(temp_dir.path(), |relative_path, _| {
        !relative_path.starts_with("input")
    })?;

    let missing_files = filter_files(
        dir_diff(
            relative(&expected_files),
            relative(&output_files_in_temp_dir),
        )?,
        |filename| !(filename.starts_with("compare_") || filename.contains("/compare_")),
    );
    let unexpected_files = dir_diff(
        relative(&output_files_in_temp_dir),
        relative(&expected_files),
    )?;
    let common_files = dir_common(
        relative(&expected_files),
        relative(&output_files_in_temp_dir),
    )?;

    let mut missmatched_files = Vec::new();

    for relative_filename in common_files {
        // Compare each file in the common filesi
        if !files_equal(
            test_case.dir.join(&relative_filename),
            temp_dir.path().join(&relative_filename),
        )
        .unwrap()
        {
            missmatched_files.push((
                temp_dir
                    .path()
                    .join(&relative_filename)
                    .to_string_lossy()
                    .to_string(),
                test_case
                    .dir
                    .join(&relative_filename)
                    .to_string_lossy()
                    .to_string(),
            ));
        }
    }
    result.missing_files = missing_files;
    result.unexpected_files = unexpected_files;
    result.mismatched_files = missmatched_files;

    // First, check all files in the temp directory that should match expected outputs

    if !(result.missing_files.is_empty()
        && result.mismatched_files.is_empty()
        && result.unexpected_files.is_empty())
    {
        // Create actual directory and copy files
        /* if actual_dir.exists() {
            fs::remove_dir_all(&actual_dir)?;
        } */
        fs::create_dir_all(&actual_dir)?;
        //copy all files from temp_dir to actual_dir
        visit_dirs(temp_dir.path(), &mut |entry| {
            let absolute_src_path = entry.path();
            let relative_src_path = absolute_src_path
                .strip_prefix(temp_dir.path())
                .context("Strip prefix from temp dir path")?;
            if absolute_src_path.is_file() {
                let dest_path = actual_dir.join(relative_src_path);
                std::fs::create_dir_all(dest_path.parent().unwrap())?;
                fs::copy(&absolute_src_path, &dest_path)?;
            }
            Ok(())
        })?;
    } else {
        //remove actual dir
        if actual_dir.exists() {
            fs::remove_dir_all(&actual_dir)?;
        }
    }
    Ok(result)
}

fn files_equal(file_a: PathBuf, file_b: PathBuf) -> Result<bool> {
    let content_a = ex::fs::read(&file_a).unwrap();
    let content_b = ex::fs::read(&file_b).unwrap();
    if content_a == content_b {
        return Ok(true);
    }
    if file_a.extension() == Some(std::ffi::OsStr::new("gz"))
        && file_b.extension() == Some(std::ffi::OsStr::new("gz"))
    {
        let uncompressed_a = read_compressed(&file_a)?;
        let uncompressed_b = read_compressed(&file_b)?;
        return Ok(uncompressed_a == uncompressed_b);
    }
    let comparison_script = file_a.with_file_name(format!(
        "compare_{}",
        file_a.file_name().unwrap().to_string_lossy()
    ));
    if comparison_script.exists() {
        let output = std::process::Command::new(&comparison_script)
            .arg(&file_a)
            .arg(&file_b)
            .output()
            .context("Failed to run comparison script. Is it executable?")?;
        if output.status.success() {
            return Ok(true);
        } else {
            let stdout = String::from_utf8_lossy(&output.stdout);
            let stderr = String::from_utf8_lossy(&output.stderr);
            println!(
                "Comparison script failed for {} with stdout: {} error: {}",
                file_a.display(),
                stdout.trim(),
                stderr.trim()
            );
            return Ok(false);
        }
    }
    Ok(false)
}

fn copy_files(input_files: &Vec<(PathBuf, String)>, target_dir: &Path) -> Result<()> {
    for (input_file, relative_path) in input_files {
        let dst_path = target_dir.join(relative_path);
        std::fs::create_dir_all(dst_path.parent().unwrap())?;
        fs::copy(input_file, &dst_path)?;
    }
    Ok(())
}

fn setup_test_environment(input_files: Vec<(PathBuf, String)>) -> Result<TempDir> {
    let temp_dir = tempfile::tempdir().context("make tempdir")?;
    copy_files(&input_files, temp_dir.path()).context("Copy input files to temp dir")?;

    Ok(temp_dir)
}

fn filter_files(files: Vec<String>, filter: impl Fn(&str) -> bool) -> Vec<String> {
    let filtered: Vec<String> = files.into_iter().filter(|f| filter(f)).collect();
    filtered
}

fn dir_diff(files_a: Vec<String>, files_b: Vec<String>) -> Result<Vec<String>> {
    let mut diff = Vec::new();
    let set_b: std::collections::HashSet<_> = files_b.iter().collect();
    for file in files_a {
        if !set_b.contains(&file) {
            diff.push(file.clone());
        }
    }
    Ok(diff)
}

fn dir_common(files_a: Vec<String>, files_b: Vec<String>) -> Result<Vec<String>> {
    let mut common = Vec::new();
    let set_b: std::collections::HashSet<_> = files_b.iter().collect();
    for file in files_a {
        if set_b.contains(&file) {
            common.push(file.clone());
        }
    }
    Ok(common)
}

fn relative(files: &[(PathBuf, String)]) -> Vec<String> {
    files.iter().map(|(_, relative)| relative.clone()).collect()
}
