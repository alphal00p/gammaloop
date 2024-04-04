#!/usr/bin/env bash

RETCODE=0;



clean_dependencies () {
    cd dependencies

    rm -f INSTALLED
    rm -f dependency_build.log
    rm -rf venv
    rm -rf symbolica
    cd fjcore
    rm -f *.o
    make clean
    cd ..

    cd ..
}

build_dependencies () {

    if test -f LOCK; then
        echo "ERROR: gammaloop dependencies are already being installed by another process. Please wait for the other process to finish before running this script again, or manually remove the dependencies/LOCK file.";
        exit 1
    fi

    echo "Run 'tail -f "$(PWD)"/dependencies/dependency_build.log' to follow installation progress";

    cd dependencies
    touch LOCK
    RETCODE=$RETCODE+$?
    if [ ! $(($RETCODE)) == 0 ]
    then
        echo "ERROR: You do not have write permissions to the dependencies directory. Please run the script with the appropriate permissions.";
        exit $(($RETCODE))
    fi
    
    rm -f dependency_build.log

    if ! test -d venv; then
        PYTHON3BIN=$(which python3)
        echo "Setting up Python venv with "$PYTHON3BIN" ...";
        $PYTHON3BIN -m venv venv --system-site-packages --prompt gammaloop_venv
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            echo "ERROR: could not create Python venv";
            rm -f LOCK
            exit $(($RETCODE))
        fi
        source venv/bin/activate
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            echo "ERROR: could not activate Python venv";
            rm -f LOCK
            exit $(($RETCODE))
        fi
        echo "Installing Python dependencies in venv ...";
        pip install --upgrade pip &>> dependency_build.log
        pip install maturin &>> dependency_build.log
        pip install -r ../requirements.txt &>> dependency_build.log
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            echo "ERROR: could not install python dependencies";
            rm -f LOCK
            exit $(($RETCODE))
        fi
        python -c 'import site; print(site.getsitepackages())' > venv/site_paths.txt
    else
        source venv/bin/activate
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            echo "ERROR: failed to activate Python venv";
            rm -f LOCK
            exit $(($RETCODE))
        fi
    fi
    
    if ! test -d symbolica; then
        echo "Cloning symbolica ...";
        ${CMD_TO_ACCESS_SYMBOLICA:-git clone https://github.com/alphal00p/symbolica}
        $CMD_TO_ACCESS_SYMBOLICA &>> dependency_build.log
    fi

    if ! test -f symbolica/symbolica_path.txt; then
        cd symbolica
        echo "Building symbolica ...";
        maturin develop --release &>> ../dependency_build.log
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            echo "ERROR: failed to install symbolica. Check the logs in dependencies/dependency_build.log for more information.";
            rm -f LOCK
            exit $(($RETCODE))
        fi
        cd ..
    fi 

    if ! test -f fjcore/libfjcore.a; then
        echo "Building fjcore ...";
        cd fjcore
        make -j8 &>> ../dependency_build.log
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            echo "ERROR: failed to install fjcore..  Check the logs in dependencies/dependency_build.log for more information.";
            rm -f LOCK
            exit $(($RETCODE))
        fi
        cd ..
    fi

    touch INSTALLED
    rm -f LOCK
    
    echo "All dependencies installed successfully.";
    cd ..
}

if [ $# -eq 0 ]
    then
        echo "Building all dependencies...";
        build_dependencies
elif [ "$1" == "clean" ]
    then
        echo "Cleaning all dependencies...";
        clean_dependencies
else
    echo "Invalid argument";
	RETCODE=1;
fi

exit $(($RETCODE))
