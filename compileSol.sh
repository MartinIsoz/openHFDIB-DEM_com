#!/bin/sh

rootDir=$(pwd)

echo "Working from directory: $rootDir"

# directories to skip during compilation
skip_dirs="VOFHFDIB multiphaseInterHFDIBFoam"

# compile the solvers
cd applications/solvers
for pd in ./*/ ; do
    [ -L "${pd%/}" ] && continue
    cd $pd
    for sd in */ ; do
        [ -L "${sd%/}" ] && continue
        dir_name=$(basename "$sd")
        for skip in $skip_dirs; do
            if [ "$dir_name" = "$skip" ]; then
                continue 2
            fi
        done
        cd $sd
        wclean
        wmake
        cd ..
    done
    cd ..
done
cd $rootDir

exit 0
