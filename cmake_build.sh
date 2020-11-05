#!/bin/bash

echo "Initialize Submodules"
git submodule init
git submodule update

echo "Installing CNPY for storing data!"
mkdir ./third-party/cnpy/build
cd ./third-party/cnpy/build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8
sudo make install

cd ../../..

echo "Installing ModernRootics for IK!"
mkdir ./third-party/ModernRobotics/build
cd ./third-party/ModernRobotics/build
cmake -DCMAKE_BUILD_TYPE=Release ..
make -j8 
sudo make install

cd ../../..
echo "Installing admm standalone"
mkdir ./standalone/build
cd ./standalone/build
cmake -DCMAKE_BUILD_TYPE=Release -DCNPY_DIR="../third-party/cnpy/cmake" -DModernRoboticsCpp_DIR="../third-party/ModernRobotics/cmake" ..
make -j8
sudo make install

cd ../../

echo "Installation Complete"



