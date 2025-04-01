#!/usr/bin/env bash
echo "Building cli binary... ["$@"]"
cargo build "$@" --bin cli --features="binary" --no-default-features;
exit $?
