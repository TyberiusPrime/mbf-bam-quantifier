
use human_panic::{setup_panic, Metadata};
use std::path::PathBuf;

use anyhow::{Context, Result};

fn print_usage(exit_code: i32) -> ! {
    let this = std::env::args().next().unwrap();
    eprintln!(
        "Usage:
    {this} <config.toml>
    {this} --version # output version and exit(0)

",
    );
    std::process::exit(exit_code);
}

fn main() -> Result<()> {
    setup_panic!(
        Metadata::new(env!("CARGO_PKG_NAME"), env!("CARGO_PKG_VERSION"))
            //.authors("My Company Support <support@mycompany.com>")
            .homepage("https://github.com/TyberiusPrime/mbf-bam-quantifier")
            .support("Open a github issue at https://github.com/TyberiusPrime/mbf-bam-quantifier/issues/new and attach the crash report.")
    );

    assert!(
        !std::env::args().any(|x| x == "--test-friendly-panic"),
        "friendly panic test!"
    );

    if std::env::args().any(|x| x == "--version") {
        println!("{}", env!("CARGO_PKG_VERSION"));
        std::process::exit(0);
    }
    if std::env::args().any(|x| x == "--help") {
        print_usage(0);
    }
    if std::env::args().len() < 2 {
        print_usage(1);
    }
    let toml_file = std::env::args()
        .nth(1)
        .context("First argument must be a toml file path.")?;
    let toml_file = PathBuf::from(toml_file);
    /* let current_dir = std::env::args()
        .nth(2)
        .map_or_else(|| std::env::current_dir().unwrap(), PathBuf::from); */
    if let Err(e) = mbf_bam_quantifier::run(&toml_file) {
        eprintln!(
            "Unfortunatly an error was detected and lead to an early exit.\n\nDetails: {e:?}",
        );
        std::process::exit(1);
    }
    Ok(())
}
