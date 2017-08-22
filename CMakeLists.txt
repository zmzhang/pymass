cmake_minimum_required(VERSION 2.8.0)

project(pymass)

# set module path for eigen, numpy
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}  "../thirdparty/eigen")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}  "../thirdparty/expat")
set(CMAKE_PREFIX_PATH ${CMAKE_PREFIX_PATH}  "../thirdparty/tbb")

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


# find base64SIMD
get_filename_component(BASE64_ABS_PATH "../thirdparty/base64SIMD" ABSOLUTE)
include_directories("${BASE64_ABS_PATH}/include")
link_directories("${BASE64_ABS_PATH}/lib")

# find tbb
find_package(tbb) 


find_package(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

add_definitions(-DUNICODE -D_UNICODE -D_CRT_SECURE_NO_DEPRECATE /wd4251 /wd4819 /wd4244)


include_directories(cpp)


add_subdirectory(cpp)
add_subdirectory(python)
