use anyhow::{bail, Context, Result};
use std::fs::{self, DirEntry};
use std::io::Read;
use std::os::unix::fs::PermissionsExt;
use std::path::{Path, PathBuf};
use std::process;
use tempfile::TempDir;

const CLI_UNDER_TEST: &str = "mbf-bam-quantifier";

fn main() -> Result<()> {
    human_panic::setup_panic!();
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
        let repeat_count = fs::read_to_string(test_case.dir.join("repeat"))
            .map(|x| {
                x.trim()
                    .parse::<usize>()
                    .expect("Repeat file with non number")
            })
            .unwrap_or(1);

        for repeat in 0..repeat_count {
            print!("\n  Running test: {} {}", test_case.dir.display(), repeat);
            let start = std::time::Instant::now();
            let test_result = if test_case.is_panic {
                run_panic_test(&test_case, processor_path.as_ref())
            } else {
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

fn perform_test(test_case: &TestCase, processor_cmd: &Path) -> Result<TestOutput> {
    let mut result = TestOutput {
        stdout: String::new(),
        stderr: String::new(),
        return_code: 0,
        missing_files: Vec::new(),
        mismatched_files: Vec::new(),
        unexpected_files: Vec::new(),
    };

    let temp_dir = setup_test_environment(&test_case.dir).context("Setup test dir")?;

    // Run the processor
    let config_file = temp_dir.path().join("input.toml");
    //chdir to temp_dir

    let actual_dir = test_case.dir.join("actual");
    // Create actual directory and copy files
    if actual_dir.exists() {
        fs::remove_dir_all(&actual_dir)?;
    }
    fs::create_dir_all(&actual_dir)?;
    //copy all files from temp_dir to actual_dir
    for entry in fs::read_dir(temp_dir.path())? {
        let entry = entry?;
        let src_path = entry.path();
        if src_path.is_file() {
            let dest_path = actual_dir.join(src_path.file_name().unwrap());
            fs::copy(&src_path, &dest_path)?;
        }
    }

    let proc = std::process::Command::new(processor_cmd)
        .arg(&config_file)
        .arg(temp_dir.path())
        .current_dir(temp_dir.path())
        .output()
        .context(format!("Failed to run {CLI_UNDER_TEST}"))?;

    let stdout = String::from_utf8_lossy(&proc.stdout);
    let stderr = String::from_utf8_lossy(&proc.stderr);
    result.return_code = proc.status.code().unwrap_or(-1);
    result.stdout = stdout.to_string();
    result.stderr = stderr.to_string();

    //for comparison
    fs::write(temp_dir.path().join("stdout"), stdout.as_bytes())
        .context("Failed to write stdout to file")?;
    /* fs::write(temp_dir.path().join("stderr"), stderr.as_bytes())
    .context("Failed to write stderr to file")?; */
    //for debugging..
    fs::write(actual_dir.as_path().join("stdout"), stdout.as_bytes())
        .context("Failed to write stdout to file")?;
    fs::write(actual_dir.as_path().join("stderr"), stderr.as_bytes())
        .context("Failed to write stderr to file")?;

    // First, check all files in the temp directory that should match expected outputs
    visit_dirs(temp_dir.path(), &mut |entry: &DirEntry| -> Result<()> {
        let path = entry.path();
        let relative_path = path
            .strip_prefix(temp_dir.path())
            .context("Strip prefix from temp dir path")?;

        // Skip input files and special files
        if let Some(file_name) = path.file_name() {
            let file_name_str = file_name.to_string_lossy();
            let parent_name = path
                .parent()
                .unwrap()
                .file_name()
                .unwrap()
                .to_string_lossy();
            if file_name_str.starts_with("input")
                || file_name_str.starts_with("ignore_")
                || parent_name.starts_with("ignore_")
                || file_name_str.starts_with("ignore_")
            {
                return Ok(());
            }
            for only_if_expected_filename in ["stdout"] {
                if file_name_str == only_if_expected_filename {
                    //only check if there's an expected stdout
                    let expected_path = test_case.dir.join(only_if_expected_filename);
                    if !expected_path.exists() {
                        // If there's no expected stdout, skip this file
                        return Ok(());
                    }
                }
            }
        }

        if path.is_file() {
            let expected_path = test_case.dir.join(relative_path);
            if expected_path.exists() {
                // Compare files
                let expected_content = fs::read(&expected_path)?;
                let actual_content = fs::read(&path)?;

                if expected_content != actual_content {
                    //if compressed, compare uncompressed
                    if expected_path
                        .extension()
                        .map_or(false, |ext| ext == "gz" || ext == "zst")
                    {
                        let expected_uncompressed = read_compressed(&expected_path)?;
                        let actual_uncompressed = read_compressed(&path)?;
                        if expected_uncompressed != actual_uncompressed {
                            result.mismatched_files.push((
                                path.to_string_lossy().to_string(),
                                expected_path.to_string_lossy().to_string(),
                            ));
                        }
                    } else if expected_path.extension().map_or(false, |ext| ext == "json") {
                        //we need to avoid the <working_dir> in reports
                        let actual_content = std::str::from_utf8(&actual_content)
                            .context("Failed to convert actual content to string")?
                            .replace(temp_dir.path().to_string_lossy().as_ref(), "WORKINGDIR")
                            .as_bytes()
                            .to_vec();
                        //support for _internal_read_count checks.
                        //thease are essentialy <=, but we just want to compare json as strings, bro
                        if actual_content != expected_content {
                            fs::write(&path, &actual_content)
                                .context("Failed to write actual content to file")?;
                            result.mismatched_files.push((
                                path.to_string_lossy().to_string(),
                                expected_path.to_string_lossy().to_string(),
                            ));
                        }
                    } else if expected_path
                        .extension()
                        .map_or(false, |ext| ext == "progress")
                    {
                        //remove all numbres from actual and expected and compare again
                        let expected_wo_numbers = regex::Regex::new(r"\d+")
                            .unwrap()
                            .replace_all(std::str::from_utf8(&expected_content).unwrap(), "");
                        let actual_wo_numbers = regex::Regex::new(r"\d+")
                            .unwrap()
                            .replace_all(std::str::from_utf8(&actual_content).unwrap(), "");
                        if expected_wo_numbers != actual_wo_numbers {
                            result.mismatched_files.push((
                                path.to_string_lossy().to_string(),
                                expected_path.to_string_lossy().to_string(),
                            ));
                        }
                    } else {
                        result.mismatched_files.push((
                            path.to_string_lossy().to_string(),
                            expected_path.to_string_lossy().to_string(),
                        ));
                    }
                }
            } else {
                // Expected file doesn't exist - this is a new output file
                result
                    .unexpected_files
                    .push(path.to_string_lossy().to_string());
            }
        }
        Ok(())
    })?;

    // Also check if there are any expected output files that weren't produced
    for entry in fs::read_dir(&test_case.dir)? {
        let entry = entry?;
        let expected_path = entry.path();

        if expected_path.is_file() {
            if let Some(file_name) = expected_path.file_name() {
                let file_name_str = file_name.to_string_lossy();

                // Skip non-output files
                if file_name_str.starts_with("input")
                    || file_name_str == "expected_panic.txt"
                    || file_name_str == "error"
                    || file_name_str == "repeat"
                    || file_name_str == "top.json"
                {
                    continue;
                }

                let actual_path = temp_dir.path().join(&file_name);
                if !actual_path.exists() {
                    // Expected output file was not produced
                    result
                        .missing_files
                        .push(expected_path.to_string_lossy().to_string());
                }
            }
        }
    }

    if !(result.missing_files.is_empty()
        && result.mismatched_files.is_empty()
        && result.unexpected_files.is_empty())
    {
        // Create actual directory and copy files
        if actual_dir.exists() {
            fs::remove_dir_all(&actual_dir)?;
        }
        fs::create_dir_all(&actual_dir)?;
        //copy all files from temp_dir to actual_dir
        visit_dirs(temp_dir.path(), &mut |entry| {
            let absolute_src_path = entry.path();
            let relative_src_path = absolute_src_path
                .strip_prefix(temp_dir.path())
                .context("Strip prefix from temp dir path")?;
            if absolute_src_path.is_file() {
                let dest_path = actual_dir.join(&relative_src_path);
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

fn setup_test_environment(test_dir: &Path) -> Result<TempDir> {
    let temp_dir = tempfile::tempdir().context("make tempdir")?;

    // Copy input.toml
    let input_toml_src = test_dir.join("input.toml");
    let input_toml_dst = temp_dir.path().join("input.toml");
    fs::copy(&input_toml_src, &input_toml_dst).context("copy input file")?;

    // Copy any input*.fq* files
    for entry in fs::read_dir(test_dir)? {
        let entry = entry?;
        let path = entry.path();
        if path.is_file() {
            if let Some(file_name) = path.file_name() {
                let file_name_str = file_name.to_string_lossy();
                if file_name_str.starts_with("input_") {
                    let dst_path = temp_dir.path().join(file_name);
                    fs::copy(&path, &dst_path)?;
                }
            }
        }
    }

    Ok(temp_dir)
}
