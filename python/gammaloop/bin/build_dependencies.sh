#!/usr/bin/env bash

RETCODE=0;



clean_dependencies () {
    cd dependencies

    rm -f INSTALLED
    rm -f LOCK
    rm -f dependency_build.log
    rm -f test_quad_math/test_quad_math
    rm -rf venv
    rm -rf symbolica
    cd fjcore
    rm -f *.o
    make clean >> /dev/null 2>&1
    cd ..

    cd ..
}

build_dependencies () {

    if test -f LOCK; then
        echo -e "\033[91mERROR: gammaloop dependencies are already being installed by another process. Please wait for the other process to finish before running this script again, or manually remove the dependencies/LOCK file.\033[0m";
        exit 1
    fi

    echo "Run 'tail -f "$(pwd)"/dependencies/dependency_build.log' to follow installation progress.";

    cd dependencies
    touch LOCK
    RETCODE=$RETCODE+$?
    if [ ! $(($RETCODE)) == 0 ]
    then
        echo -e "\033[91mERROR: You do not have write permissions to the dependencies directory. Please run the script with the appropriate permissions.\033[0m";
        exit $(($RETCODE))
    fi
    
    rm -f dependency_build.log

    FINALMESSAGE="\033[92mAll dependencies installed successfully.\033[0m";

    CPPCOMPILER="${CXX:-g++}"
    CCOMPILER="${CC:-cc}"
    # We must also test explictly cc as it is used *explicitely* when building some of quadruple precision rust crates
    FORCEDCCOMPILER="cc"
#    if ! test -f test_quad_math/test_quad_math; then
    if ! true; then
        cd test_quad_math

        echo "Testing quadruple precision support with C++ compiler "$CPPCOMPILER" ...";
        $CPPCOMPILER test_quad_math.c -o test_quad_math -lquadmath >> ../dependency_build.log 2>&1
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            cat ../dependency_build.log
            echo -e "\033[91mERROR: could not compile with quadruple precision with compiler "$CPPCOMPILER". Make sure you are using GNU GCC and not clang.\033[0m"
            rm -f LOCK
            exit $(($RETCODE))
        fi
        ./test_quad_math >> ../dependency_build.log 2>&1
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            cat ../dependency_build.log
            echo -e "\033[91mERROR: could not run code testing quadruple precision compiled with compiler "$CPPCOMPILER".\033[0m"
            rm -f LOCK
            exit $(($RETCODE))
        fi
        rm -f test_quad_math

        echo "Testing quadruple precision support with C++ compiler "$CCOMPILER" ...";
        $CCOMPILER test_quad_math.c -o test_quad_math -lquadmath >> ../dependency_build.log 2>&1
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            cat ../dependency_build.log
            echo -e "\033[91mERROR: could not compile with quadruple precision with compiler "$CCOMPILER". Make sure you are using GNU GCC and not clang.\033[0m"
            rm -f LOCK
            exit $(($RETCODE))
        fi
        ./test_quad_math >> ../dependency_build.log 2>&1
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            cat ../dependency_build.log
            echo -e "\033[91mERROR: could not run code testing quadruple precision compiled with compiler "$CCOMPILER".\033[0m"
            rm -f LOCK
            exit $(($RETCODE))
        fi

        if [ ! $FORCEDCCOMPILER == $CCOMPILER ]
        then
            rm -f test_quad_math

            echo "Testing quadruple precision support with C++ compiler "$FORCEDCCOMPILER" ...";
            $FORCEDCCOMPILER test_quad_math.c -o test_quad_math -lquadmath >> ../dependency_build.log 2>&1
            RETCODE=$RETCODE+$?
            if [ ! $(($RETCODE)) == 0 ]
            then
                cat ../dependency_build.log
                echo -e "\033[91mERROR: could not compile with quadruple precision with compiler "$FORCEDCCOMPILER". Make sure you are using GNU GCC and not clang.\033[0m"
                rm -f LOCK
                exit $(($RETCODE))
            fi
            ./test_quad_math >> ../dependency_build.log 2>&1
            RETCODE=$RETCODE+$?
            if [ ! $(($RETCODE)) == 0 ]
            then
                cat ../dependency_build.log
                echo -e "\033[91mERROR: could not run code testing quadruple precision compiled with compiler "$FORCEDCCOMPILER".\033[0m"
                rm -f LOCK
                exit $(($RETCODE))
            fi
        fi
        
        cd ..
    fi

    PYTHON3BIN=$(which python3)
    if [ "$1" == "with_venv" ]
    then

        if ! test -d venv; then
            echo "Setting up Python venv with "$PYTHON3BIN" ...";
            $PYTHON3BIN -m venv venv --system-site-packages --prompt gammaloop_venv
            RETCODE=$RETCODE+$?
            if [ ! $(($RETCODE)) == 0 ]
            then
                echo -e "\033[91mERROR: could not create Python venv\033[0m";
                rm -f LOCK
                exit $(($RETCODE))
            fi
            source venv/bin/activate
            PYTHON3BIN="python"
            RETCODE=$RETCODE+$?
            if [ ! $(($RETCODE)) == 0 ]
            then
                echo -e "\033[91mERROR: could not activate Python venv\033[0m";
                rm -f LOCK
                exit $(($RETCODE))
            fi
            echo "Installing Python dependencies in venv ...";
            pip install --upgrade pip >> dependency_build.log 2>&1
            RETCODE=$RETCODE+$?
            pip install maturin >> dependency_build.log 2>&1
            RETCODE=$RETCODE+$?
            pip install -r ../requirements.txt >> dependency_build.log 2>&1
            RETCODE=$RETCODE+$?
            if [ ! $(($RETCODE)) == 0 ]
            then
                echo -e "\033[91mERROR: could not install python dependencies\033[0m";
                rm -f LOCK
                exit $(($RETCODE))
            fi
            $PYTHON3BIN -c 'import site; print(site.getsitepackages())' > venv/site_paths.txt
        else
            echo "Activating Python venv with "$(PWD)"/venv/bin/activate ...";
            source venv/bin/activate
            RETCODE=$RETCODE+$?
            PYTHON3BIN="python"
            if [ ! $(($RETCODE)) == 0 ]
            then
                echo -e "\033[91mERROR: failed to activate Python venv\033[0m";
                rm -f LOCK
                exit $(($RETCODE))
            fi
        fi

    else
        $PYTHON3BIN check_python_dependencies.py >> dependency_build.log 2>&1
        if [ ! $(($?)) == 0 ]
        then
            echo "Installing Python dependencies with pip...";
            $PYTHON3BIN -m pip install -r ../requirements.txt >> dependency_build.log 2>&1
            if [ ! $(($?)) == 0 ]
            then
                FINALMESSAGE="\033[92mAll non Python dependencies installed successfully.\033[0m";
                FINALMESSAGE=$FINALMESSAGE"\n\033[93mWARNING: could not install Python dependencies with pip. You will need to install them manually with '"$PYTHON3BIN" -m pip install -r requirements.txt'. You can also consider doing this within a virtual environment with '"$PYTHON3BIN" -m venv .venv'\033[0m"
                echo -e "\033[93mWARNING: could not install Python dependencies with pip. You will need to install them manually with '"$PYTHON3BIN" -m pip install -r requirements.txt'. You can also consider doing this within a virtual environment with '"$PYTHON3BIN" -m venv .venv'\033[0m";
            fi
        else
            echo "All Python dependencies already installed.";
        fi
    fi

    if ! test -d symbolica; then
        CMD_TO_ACCESS_SYMBOLICA="${CMD_TO_ACCESS_SYMBOLICA:-git clone -b symbolica_fork_for_v0_3_3 https://github.com/alphal00p/symbolica}"
        echo "Cloning symbolica with '"$CMD_TO_ACCESS_SYMBOLICA"' ...";
        $CMD_TO_ACCESS_SYMBOLICA >> dependency_build.log 2>&1
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            cat ../dependency_build.log;
            echo -e "\033[91mERROR: failed to clone symbolica with command '"$CMD_TO_ACCESS_SYMBOLICA"'. Check the logs in dependencies/dependency_build.log for more information.\033[0m";
            rm -f LOCK
            exit $(($RETCODE))
        fi
        echo "Patching symbolica ...";
        cd symbolica && ln -s ../../bin/patch_symbolica.py . && ./patch_symbolica.py >> ../dependency_build.log 2>&1
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            cat ../dependency_build.log;
            echo -e "\033[91mERROR: failed to patch symbolica. Check the logs in dependencies/dependency_build.log for more information.\033[0m";
            rm -f LOCK
            exit $(($RETCODE))
        fi
        cd ..
    fi
    
	if ! test -f symbolica/symbolica_path.txt; then
        cd symbolica
        
		SYMBOLICA_BUILD_PROFILE="${SYMBOLICA_BUILD_PROFILE:-release}"
        if [ "$1" == "with_venv" ]
        then

            echo "Building symbolica with maturin within a venv ...";
            maturin develop >> ../dependency_build.log 2>&1
            RETCODE=$RETCODE+$?
            if [ ! $(($RETCODE)) == 0 ]
            then
                cat ../dependency_build.log;
                echo -e "\033[91mERROR: failed to install symbolica. Check the logs in dependencies/dependency_build.log for more information.\033[0m";
                rm -f LOCK
                exit $(($RETCODE))
            fi
            $PYTHON3BIN -c "import os; import symbolica; print(os.path.abspath(os.path.join(os.path.dirname(symbolica.__file__),os.path.pardir)))" > symbolica_path.txt
            RETCODE=$RETCODE+$?
            if [ ! $(($RETCODE)) == 0 ]
            then
                cat ../dependency_build.log;
                echo -e "\033[91mERROR: failed to load symbolica Python module built by maturin. Check the logs in dependencies/dependency_build.log for more information.\033[0m";
                rm -f LOCK
                exit $(($RETCODE))
            fi

        else
            echo "Building symbolica ...";
            PYO3_PYTHON=$PYTHON3BIN cargo rustc --profile=$SYMBOLICA_BUILD_PROFILE --features=python_api --crate-type=cdylib >> ../dependency_build.log 2>&1
            RETCODE=$RETCODE+$?
            if [ ! $(($RETCODE)) == 0 ]
            then
                cat ../dependency_build.log;
                echo -e "\033[91mERROR: failed to manually build symbolica's python module. Check the logs in dependencies/dependency_build.log for more information.\033[0m";
                rm -f LOCK
                exit $(($RETCODE))
            fi
            
            if [[ $SYMBOLICA_BUILD_PROFILE = "dev" ]]
            then
                COMPILE_PATH="debug"
            else
                COMPILE_PATH=$SYMBOLICA_BUILD_PROFILE
            fi

            if test -f target/$COMPILE_PATH/libsymbolica.so; then
                ln -s target/$COMPILE_PATH/libsymbolica.so symbolica.so
            elif test -f target/$COMPILE_PATH/libsymbolica.dylib; then
                ln -s target/$COMPILE_PATH/libsymbolica.dylib symbolica.so
            elif test -f target/$COMPILE_PATH/libsymbolica.dll; then
                ln -s target/$COMPILE_PATH/libsymbolica.dll symbolica.pyd
            else
                echo -e "\033[91mERROR: failed to find manually compiled symbolica's python module. Check the logs in dependencies/dependency_build.log for more information.\033[0m";
                rm -f LOCK
                exit 1
            fi
            echo "." > symbolica_path.txt

        fi

        cd ..
    fi 

    if ! test -f fjcore/libfjcore.a; then
        echo "Building fjcore ...";
        cd fjcore
        make -j8 >> ../dependency_build.log 2>&1
        RETCODE=$RETCODE+$?
        if [ ! $(($RETCODE)) == 0 ]
        then
            cat ../dependency_build.log;
            echo -e "\033[91mERROR: failed to install fjcore. Check the logs in dependencies/dependency_build.log for more information.\033[0m";
            rm -f LOCK
            exit $(($RETCODE))
        fi
        cd ..
    fi

    touch INSTALLED
    rm -f LOCK

    echo ""
    echo "---- SUMMARY ----"
    echo ""
    
    echo -e $FINALMESSAGE;

    cd ..
}

if [ $# -eq 0 ]
    then
        echo "Building all dependencies...";
        build_dependencies no_venv
elif [ "$1" == "with_venv" ]
    then
        echo "Building all dependencies within a Python venv...";
        build_dependencies with_venv
elif [ "$1" == "clean" ]
    then
        echo "Cleaning all dependencies...";
        clean_dependencies
else
    echo "Invalid argument";
	RETCODE=1;
fi

exit $(($RETCODE))
