# Configure and build our SWIG-based python wrapper.
# We defined a helper function in swigmake/cmake/SwigPython.cmake

add_definitions(/wd4267 /wd4244)

include(SwigPython)
add_swig_python_module(TARGET pymass
                       CPLUSPLUS
                       SWIG_I_FILE ${CMAKE_CURRENT_SOURCE_DIR}/pymass.i
                       INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../cpp
                       LINK_LIBRARIES pymass
                       SWIG_INCLUDE_DIRS ${CMAKE_CURRENT_SOURCE_DIR}/../swig
                       DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/)
# If we wanted to, we could copy all the python source files to an install directory like /usr/local/lib/python2.7/site-packages
# Instead, we're just going to keep the python files where they are and install the module that swig generates right next to them. 

