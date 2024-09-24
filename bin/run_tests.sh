#!/usr/bin/env bash

RETCODE=0;

run_rust_tests() {
    export SYMBOLICA_LICENSE="${SYMBOLICA_LICENSE:-GAMMALOOP_USER}";
    cargo test --features=binary --release --no-default-features -- --test-threads=1 --skip tests_from_pytest "$@";
    # Run the tests using specific categories and showing live println()'s
    #cargo test --release --features=binary --no-default-features -- tests_inspect --test-threads=1 --nocapture
    #cargo test --release --features=binary --no-default-features -- tests_integral --test-threads=1 --nocapture
    RETCODE=$RETCODE+$?;
}

run_python_tests() {
    cd python/gammaloop;
    ./run_python_tests.sh "$@";
    RETCODE=$RETCODE+$?;
    # Pass the following option to see output on fails
    #./run_python_tests.sh "$@" -s --log-cli-level DEBUG;
    # or:
    #./run_python_tests.sh "$@" -rx;
    # Pass the following option to see output on all tests
    #./run_python_tests.sh "$@" -rP;
    cd - >& /dev/null;
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
        run_rust_tests "${@:2}"
elif [ "$1" == "python" ]
    then
        echo "Running Python tests"
        run_python_tests "${@:2}"
elif [ "$1" == "all" ]
    then
        echo "Running all tests"
        run_rust_tests "${@:2}"
        run_python_tests "${@:2}"
else
    echo "Invalid argument"
    RETCODE=1;
fi

exit $(($RETCODE))
