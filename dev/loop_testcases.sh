#/usr/bin/bash

fd "panic|\.rs\$|\.toml\$|sha256" | grep -v actual | entr cargo test --release test_case
