#!/usr/bin/env bash

RETCODE=0;

clean_dependencies () {
    cd dependencies/fjcore
    rm -f *.o
    make clean
    cd ../..
}

build_dependencies () {
    cd dependencies/fjcore
    make -j8 &> /dev/null
    RETCODE=$RETCODE+$?
    cd ../..
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
