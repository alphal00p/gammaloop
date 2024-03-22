#!/usr/bin/env bash
RETCODE=0;
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
RETCODE=$RETCODE+$?;
rm -rf $TMPDIR/test_gammaloop_deployment
mkdir $TMPDIR/test_gammaloop_deployment
mkdir $TMPDIR/test_gammaloop_deployment/wheel
maturin build --release --features "extension-module" -o /tmp/test_gammaloop_deployment/wheel
RETCODE=$RETCODE+$?;
cd $TMPDIR/test_gammaloop_deployment
python -m venv venv
source ./venv/bin/activate
python -m pip install --upgrade pip
python -m pip install pytest
python -m pip install `ls -1 $TMPDIR/test_gammaloop_deployment/wheel/*.whl`
RETCODE=$RETCODE+$?;
gammaloop --build_dependencies
RETCODE=$RETCODE+$?;
cd `ls -d1 ./venv/lib/python*/site-packages/gammaloop`
source `gammaloop -venv`
python -m pytest
RETCODE=$RETCODE+$?;
exit $(($RETCODE))
