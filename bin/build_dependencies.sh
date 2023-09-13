#!/usr/bin/env bash

clean_dependencies () {
    cd dependencies/fjcore
    rm *.o
    make clean
    cd ../..
}

build_dependencies () {
    cd dependencies/fjcore
    make -j8 >& /dev/null
    cd ../..
}

if [ $# -eq 0 ]
    then
        echo "Building all dependencies..."
        build_dependencies
elif [ "$1" == "clean" ]
    then
        echo "Cleaning all dependencies..."
        clean_dependencies
else
    echo "Invalid argument"
fi