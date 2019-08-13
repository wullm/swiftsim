#!/bin/bash

#========================================
# Script to quickly compile swift with
# the right flags and recompile as little
# as possible and necessary.
#
# use ./myinstall to just recompile using
# the previously used flags, provided the
# file .last_compile is present in this
# directory. Otherwise, it will just use
# the default flags below.
# You can also give it specific flags for
# which to compile. See below for
# possibilities (search keyword MYFLAGS)
#========================================

# change silent to false if you want the full ./configure and make outputs written
# to stdout and stderr. Otherwise, will be written to $logfile.
# script will exit if there was an error, so you'll know if something goes wrong.

silent=false
logfile=log_of_my_install
if [ -f $logfile ]; then rm $logfile; fi


DEBUGFLAGS=''
DEBUGFLAGS_IF_IN_USE="--enable-debug --enable-sanitizer --enable-optimization=no --enable-undefined-sanitizer" # if debug is selected, these debugging flags will be used.
DEFAULTFLAGS='--enable-mpi=no --disable-doxygen-doc'
DIMFLAGS='' # default 3D
GIZMOFLAGS='--with-hydro=gizmo-mfv --with-riemann-solver=hllc'
LIBFLAGS="--with-parmetis --with-jemalloc --with-hdf5=$HDF5_ROOT/bin/h5pcc"

EXTRA_CFLAGS=""


debug_program_suffix=''


function errexit() {
    # usage: errexit $? "optional message string"
    if [[ "$1" -ne 0 ]]; then
        echo "ERROR OCCURED. ERROR CODE $1"
        if [[ $# > 1 ]]; then
            echo "$2"
        fi
        exit $1
    else
        return 0
    fi
}

function file_separator() {
    # usage: filename "text to add"
    echo "=======================================================" >> $1
    echo "=======================================================" >> $1
    echo "========" $2 >> $1
    echo "=======================================================" >> $1
    echo "=======================================================" >> $1
    return 0
}




if [ ! -f ./configure ]; then
    ./autogen.sh
fi






#--------------------------------------
# standardize comp flag
#--------------------------------------

# HERE ARE MYFLAGS
case $1 in

    default | -d | d | 3 | 3d | --3d | -3d)
        echo COMPILING DEFAULT
        comp=default
    ;;

    clean | c | -c | --c | --clean)
        echo COMPILING CLEAN
        comp=clean
    ;;

    1 | --1d | -1d | -1 | --1 | 1d)
        echo COMPILING 1D SWIFT
        DIMFLAGS='--with-hydro-dimension=1'
        comp=1d
    ;;

    2 | --2d | -2d | -2 | --2 | 2d)
        echo COMPILING 2D SWIFT
        DIMFLAGS='--with-hydro-dimension=2'
        comp=2d
    ;;

    debug | deb )
        DEBUGFLAGS=$DEBUGFLAGS_IF_IN_USE
    ;;

    last | --last | -l | --l)
        echo "COMPILING WITH SAME FLAGS AS LAST TIME"
        comp=last
    ;;
    
    *)
        echo "COMPILING WITH SAME FLAGS AS LAST TIME BY WILDCARD"
        comp=last
    ;;

esac



case $2 in

    me | my | mine | debug | deb | test )
        echo "ADDING DEBUG FLAGS"
        debug_program_suffix='-debug'
        DEBUGFLAGS=$DEBUGFLAGS_IF_IN_USE
    ;;

esac


 
allflags="$LIBFLAGS ""$GIZMOFLAGS ""$DEFAULTFLAGS"" $DEBUGFLAGS"" $DIMFLAGS"" CFLAGS=$EXTRA_CFLAGS"






#--------------------------------------
# Check if reconfiguration is necessary
#--------------------------------------



reconfigure=false


if [ -f .last_compile ]; then
    last_comp=`head -n 1 .last_compile` # only first line!
    if [ $comp != "last" ]; then # if you don't just want to repeat the same compilation

        # check that you have identical flags
        # first if all in last_comp are present in allflags
        for flag in $last_comp; do
            if [[ "$allflags" != *"$flag"* ]]; then
                echo found unmatched flag $flag
                echo will reconfigure.
                reconfigure=true
                break
            fi
        done
        if [ "$reconfigure" != 'true' ]; then
            # now if all in allflags are present in last_comp
            for flag in $allflags; do
                if [[ "$last_comp" != *"$flag"* ]]; then
                    echo found unmatched flag $flag
                    echo will reconfigure.
                    reconfigure=true
                    break
                fi
            done
        fi
    else
        # if .last_compilation exists and comp is same as last, also read in the names
        allflags=$last_comp
        lastname=`sed -n 2p .last_compile`
        lastname_mpi=`sed -n 3p .last_compile`
    fi

else
    # if no .last_compile is present
    reconfigure=true
    if [ "$comp" = 'last' ]; then
        lastname=swift
        lastname_mpi=swift_mpi
    fi
fi





#-------------------------------
# configure depending on case
#-------------------------------

if [ $reconfigure = true ]; then
    file_separator $logfile "make clean"
    if [ "$silent" = 'true' ]; then
        echo make clean
        make clean >> "$logfile"
        errexit $?
    else
        make clean | tee -a $logfile
        errexit $?
    fi

    echo configure flags are:
    for flag in $allflags; do
        echo "   " $flag
    done

    ./configure $allflags

else
    echo skipping configure.
fi




#-------------------------------
# compile
#-------------------------------

if [ $silent = "true" ]; then
    echo making.
    file_separator $logfile "make"
    make -j >> $logfile
    errexit $?
else
    make -j | tee -a $logfile
    errexit $?
fi



#--------------------------------------
# store what this compilation was
#--------------------------------------
echo "$allflags" > .last_compile






#-------------------------------
# rename executables
#-------------------------------

echo renaming.

case $comp in

    default)
        execname=./examples/swift-3d"$debug_program_suffix"
        execname_mpi=./examples/swift_mpi-3d"$debug_program_suffix"
    ;;

    1d)
        execname=./examples/swift-1d"$debug_program_suffix"
        execname_mpi=./examples/swift_mpi-1d"$debug_program_suffix"
    ;;

    2d)
        execname=./examples/swift-2d"$debug_program_suffix"
        execname_mpi=./examples/swift_mpi-2d"$debug_program_suffix"
    ;;

    clean)
        execname=./examples/swift"$debug_program_suffix"
        execname_mpi=./examples/swift_mpi"$debug_program_suffix"
    ;;

    last)
        execname=$lastname
        execname_mpi=$lastname_mpi
    ;;

esac


mv ./examples/swift "$execname"
echo "./examples/swift -> $execname"
echo finished $execname
if [ -f ./examples/swift_mpi ]; then 
    mv ./examples/swift_mpi "$execname_mpi"
    echo "./examples/swift_mpi -> $execname_mpi"
    echo finished $execname_mpi
fi

# store last used names
echo "$execname" >> .last_compile
echo "$execname_mpi" >> .last_compile
