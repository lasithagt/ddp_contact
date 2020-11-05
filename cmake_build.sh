#!/bin/bash

echo "Initialize Submodules"
git submodule init
git submodule update

echo "Installing CNPY for storing data!"
mkdir ./third-party/cnpy/build
cd ./third-party/cnpy/build
cmake ..
make -j8
sudo make install

echo "Installing ModernRootics for IK!"
mkdir ./third-party/ModerRobotics/build
cd ./third-party/ModernRobotics/build
cmake ..
make -j8
sudo make install

echo "Installing admm standalone"
mkdir ./standalone/build
cd ./standalone/build
cmake ..
make -j8
sudo make install

echo "Installation Complete"



