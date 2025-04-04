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
echo "Local build for maturin wheel"
./bin/build_dependencies.sh
RETCODE=$RETCODE+$?;
rm -rf $TMPDIR/test_gammaloop_deployment
mkdir $TMPDIR/test_gammaloop_deployment
mkdir $TMPDIR/test_gammaloop_deployment/wheel
maturin build --release --features "extension-module" -o $TMPDIR/test_gammaloop_deployment/wheel
RETCODE=$RETCODE+$?;
cd $TMPDIR/test_gammaloop_deployment
echo "Creating virtual enviroment for testing deployment"
virtualenv venv
source ./venv/bin/activate
python3 -m pip install --upgrade pip
python3 -m pip install pytest
python3 -m pip install `ls -1 $TMPDIR/test_gammaloop_deployment/wheel/*.whl`
RETCODE=$RETCODE+$?;
echo "Building dependencies in test deployment"
gammaloop --build_dependencies
RETCODE=$RETCODE+$?;
cd `ls -d1 ./venv/lib/python*/site-packages/gammaloop`
RETCODE=$RETCODE+$?;
#source `gammaloop -venv`
echo "Running gammaloop tests withing deployed environment"
python3 -m pytest --max-runtime 15.0 -k "test_scalar_triangle"
RETCODE=$RETCODE+$?;
if [ $(($RETCODE)) == 0 ]
then
    echo -e "\033[92mDEPLOYMENT TEST SUCCESSFUL\033[0m";	
else
    echo -e "\033[91mDEPLOYMENT TEST FAILED\033[0m";
fi
exit $(($RETCODE))
