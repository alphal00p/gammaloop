#!/usr/bin/env bash
echo "Building cli binary..."
cargo build --release "$@" --bin cli --features="binary" --no-default-features;
exit $?
