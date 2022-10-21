@echo off
pushd %~dp0\Surface_Framework_Cmake\
mkdir build
cd ./build
cmake ..
echo [Done] Build finish. Please refer to the directory ./Surface_Framework_Camke/build/SurfaceFrameworkCmake.sln
cd ..
popd