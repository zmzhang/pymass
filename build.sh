mkdir thirdparty && cd thirdparty

git@github.com:RLovelett/eigen.git && cd eigen
git checkout tags/3.3.3 && cd ..

wget https://github.com/01org/tbb/releases/download/2017_U6/tbb2017_20170412oss_lin.tgz
tar xzvf tbb2017_20170412oss_lin.tgz
mv tbb2017_20170412oss tbb

git clone git@github.com:BurningEnlightenment/base64-cmake.git && cd base64-cmake
mkdir build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=../../base64SIMD && make
make install
cd .. && cd ..

git clone https://github.com/libexpat/libexpat.git
cd libexpat
git checkout tags/R_2_2_4
wget https://raw.githubusercontent.com/zmzhang/pymass/master/patches/libexpat.patch
git apply libexpat.patch
mkdir build && cd build
cmake ../expat -DCMAKE_INSTALL_PREFIX=../../expat && make && make install
cd .. && cd .. && cd ..


sudo apt-get update
sudo apt-get upgrade
sudo apt install python-dev
sudo apt install python-pip
pip install numpy


git clone git@github.com:zmzhang/pymass.git && cd pymass
mkdir build && cd build
cmake .. && make






