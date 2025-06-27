#/usr/bin/bash

fd "panic|\.rs\$|\.toml\$|sha256" | grep -v actual | entr cargo run --release --bin mbf-bam-quantifier-test-runner $@
