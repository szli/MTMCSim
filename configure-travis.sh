mkdir build_travis
cd build_travis
cmake -DCMAKE_BUILD_TYPE=Debug ../ 2>&1 | tee configure.out