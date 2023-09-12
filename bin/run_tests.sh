#!/usr/bin/env bash

run_rust_tests () {
    cargo test --release --features=binary --no-default-features -- --test-threads=1
    # Run the tests using specific categories and showing live println()'s
    #cargo test --release --features=binary --no-default-features -- tests_inspect --test-threads=1 --nocapture
    #cargo test --release --features=binary --no-default-features -- tests_integral --test-threads=1 --nocapture   
}
run_python_tests () {
    cd python/gammaloop;
    ./run_python_tests.sh;
    cd -;
}

if [ $# -eq 0 ]
    then
        echo "Running all tests"
        run_rust_tests
        run_python_tests
elif [ "$1" == "rust" ]
    then
        echo $1
        echo "Running Rust tests"
        run_rust_tests
elif [ "$1" == "python" ]
    then
        echo "Running Python tests"
        run_python_tests
elif [ "$1" == "all" ]
    then
        echo "Running all tests"
        run_rust_tests
        run_python_tests
else
    echo "Invalid argument"
fi
