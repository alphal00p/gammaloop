#!/usr/bin/env bash
cargo build $@ --bin cli --features="binary" --no-default-features;