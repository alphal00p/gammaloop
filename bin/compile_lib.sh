#!/usr/bin/env bash
RETCODE=0;

# Check if --stable-abi flag is passed
STABLE_ABI=false
CARGO_ARGS=()
for arg in "$@"; do
    case $arg in
        --stable-abi)
            STABLE_ABI=true
            ;;
        *)
            CARGO_ARGS+=("$arg")
            ;;
    esac
done

if [ "$STABLE_ABI" = true ]; then
    echo "Building pyo3 library with stable ABI... [${CARGO_ARGS[@]}]"
    cargo build "${CARGO_ARGS[@]}" --lib --features="python_abi";
else
    echo "Building pyo3 library... [${CARGO_ARGS[@]}]"
    cargo build "${CARGO_ARGS[@]}" --lib --features="python_api";
fi
RETCODE=$RETCODE+$?;
rm -f ./python/gammaloop/_gammaloop.so;
if [[ " ${CARGO_ARGS[@]} " =~ " --release " ]] || [[ " ${CARGO_ARGS[@]} " =~ " --profile=release " ]]
then
    COMPILE_PATH="release"
else
    COMPILE_PATH="debug"
fi
if [ "$(uname)" == "Darwin" ]; then
    ln -s ../../target/$COMPILE_PATH/lib_gammaloop.dylib ./python/gammaloop/_gammaloop.so;
	RETCODE=$RETCODE+$?;
else
    ln -s ../../target/$COMPILE_PATH/lib_gammaloop.so ./python/gammaloop/_gammaloop.so;
	RETCODE=$RETCODE+$?;
fi
if [ ! -f ./python/gammaloop/_gammaloop.so ]; then
  echo -e "\n\033[0;31mERROR: Could not find target pyo3 library in rust build '$COMPILE_PATH'.\033[0m\n"
  RETCODE=$RETCODE+1
fi

exit $(($RETCODE))
