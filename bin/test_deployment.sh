#!/usr/bin/env bash
TMPDIR="/tmp"
echo "Testing python wheel pip deployment in "$TMPDIR"/test_gammaloop_deployment..."
echo "The gammaloop module will be installed by pip in "$TMPDIR"/test_gammaloop_deployment/venv"
if [ "$1" == "clean" ]
    then
        echo "Cleaning project"
        cargo clean
        ./bin/build_dependencies.sh clean
fi
./bin/build_dependencies.sh
rm -rf $TMPDIR/test_gammaloop_deployment
mkdir $TMPDIR/test_gammaloop_deployment
mkdir $TMPDIR/test_gammaloop_deployment/wheel
maturin build --release -o /tmp/test_gammaloop_deployment/wheel
cd $TMPDIR/test_gammaloop_deployment
python -m venv venv
./venv/bin/python -m pip install --upgrade pip
./venv/bin/python -m pip install pytest
./venv/bin/python -m pip install `ls -1 $TMPDIR/test_gammaloop_deployment/wheel/*.whl`
cd `ls -d1 ./venv/lib/python*/site-packages/gammaloop`
../../../../bin/python -m pytest