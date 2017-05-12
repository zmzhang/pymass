cmake_minimum_required(VERSION 2.8.0)

project(pymass)

# set module path for eigen, numpy
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}  "../thirdparty/eigen")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}  "../thirdparty/expat")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}  "../thirdparty/boost")

# find eigen library
find_package( Eigen3 )
if (EIGEN3_FOUND)
  include_directories(${EIGEN3_INCLUDE_DIR} )
else()
  message(FATAL_ERROR "Could not find eigen, which is a required depedency")
endif()


# find EXPAT library
find_package(EXPAT REQUIRED)
if (EXPAT_FOUND)
  include_directories(${EXPAT_INCLUDE_DIRS})
else()
  message(FATAL_ERROR "Could not find expat, which is a required depedency")
endif()


# find boost
set(Boost_USE_STATIC_LIBS OFF) 
set(Boost_USE_MULTITHREADED ON)  
set(Boost_USE_STATIC_RUNTIME OFF) 
find_package(Boost) 
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS}) 
else()
  message(FATAL_ERROR "Could not find BOOST, which is a required depedency")
endif()

# find base64SIMD
get_filename_component(BASE64_ABS_PATH "../thirdparty/base64SIMD" ABSOLUTE)
include_directories("${BASE64_ABS_PATH}/include")
link_directories("${BASE64_ABS_PATH}/lib")

# find CUDA
FIND_PACKAGE(CUDA)
if(CUDA_FOUND)
	message("CUDA is found.")
	set(CUDA_NVCC_FLAGS "-Xcompiler \"/wd 4819 /wd 4244 /wd 4267\" ")
	set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; -O3 -gencode arch=compute_52,code=sm_52)
	set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}; --expt-relaxed-constexpr)
else(CUDA_FOUND)
    message("CUDA is not installed on this system.")
endif()


add_definitions(-DUNICODE -D_UNICODE -D_CRT_SECURE_NO_DEPRECATE /wd4251 /wd4819 /wd4244)


include_directories(cpp)
include_directories(cuda)


add_subdirectory(cpp)
add_subdirectory(python)
add_subdirectory(cuda)
