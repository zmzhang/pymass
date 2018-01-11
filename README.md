# PyMass
Package for analyzing MS with Python

It can provide the following functionalities now:


* mzXMLParser for fast and efficient mzXML parse
* FPIC method for extracting PICs from raw LC-MS dataset effectively and quickly


In future, more file formats will be supported and more methods will be implemented into PyMass package, so researchers can create complex analysis workflows for LC-MS datasets in Python with ease.

# Install

## Required Dependencies

* Hardware
	* [Modern CPUs](https://en.wikipedia.org/wiki/Advanced_Vector_Extensions#CPUs_with_AVX2) with Advanced Vector Extensions 2(AVX2) instructions 

* Windows:
	* [Visual Studio Community 2015 with Update 3](http://download.microsoft.com/download/b/e/d/bedddfc4-55f4-4748-90a8-ffe38a40e89f/vs2015.3.com_enu.iso)
	* [Anaconda Python 3.6.0 64bit](https://repo.continuum.io/archive/Anaconda3-4.3.1-Windows-x86_64.exe)
	* [SWIG 3.0.10](https://sourceforge.net/projects/swig/files/swigwin/swigwin-3.0.10/)
	* [CMake 3.7.1](https://cmake.org/files/v3.7/cmake-3.7.1-win64-x64.msi)
	* [Eigen 3.3.3](http://bitbucket.org/eigen/eigen/get/3.3.3.zip) 
	* [Threading Building Blocks 2017 Update 6](https://github.com/01org/tbb/releases/download/2017_U6/tbb2017_20170412oss_win.zip)
	* [libexpat 2.2.4](https://github.com/libexpat/libexpat/archive/R_2_2_4.tar.gz)
	* [base64SIMD](https://github.com/BurningEnlightenment/base64-cmake)
	
* Linux:
	* Ubuntu 16.04
	* GCC
	* Python
	* CMake
	* SWIG

## Download

* Download [pymass](https://github.com/zmzhang/pymass/archive/master.zip)
* Unzip it into pymass directory

## Compile
* Windows:
	* Open "VS2015 x64 Native Tools Command Prompt" 
	* Run following commands in the prompt

		```shell
		cd pymass
		mkdir build
		cd build
		cmake .. -G "NMake Makefiles" -DCMAKE_BUILD_TYPE=Release
		nmake
		nmake install
		```
* Linux:
	* PyMass can be built and run smoothly in Ubuntu Linux 16.04. We provide a bash script to download thirdparty libraries, apply the patches, build pymass automatically
		```shell
		wget https://github.com/zmzhang/pymass/raw/master/build.sh
		chmod +x build.sh
		./build.sh
		```	
# Usage

* Go to pymass/python directory
* Download MM14 dataset from this [url](https://msbi.ipb-halle.de/download/Sample-1.tar.bz2) and unzip it
* Run following Python code fragment to parse mzXML file and extract PICs from it

	```python
	from _pymass import mzXMLParser, FPICs
	import sys
	mzfile="MM14_20um.mzxml"
	mzfile=mzfile.encode(sys.getfilesystemencoding())
	parser=mzXMLParser()
	lcms = parser.parseFile(mzfile)
	pics = FPICs(lcms, 300.0, 100.0, 0.5)
	```

# Contact

For any questions, please contact:

[zmzhang@csu.edu.cn](mailto:zmzhang@csu.edu.cn)
