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
  include_directories( ${EIGEN3_INCLUDE_DIR} )
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



add_definitions(-DUNICODE -D_UNICODE -D_CRT_SECURE_NO_DEPRECATE)


add_subdirectory(cpp)
add_subdirectory(python)
#target_link_libraries(${PROJECT_NAME} ${EXPAT_LIBRARIES})
