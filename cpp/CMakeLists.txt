include(GenerateExportHeader)
add_library(pymass SHARED utils.cpp LCMS.cpp mzXMLParser.cpp)
target_include_directories(pymass PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries(pymass ${EXPAT_LIBRARIES} base64
								debug     ${TBB_LIBRARIES_DEBUG}
								optimized ${TBB_LIBRARIES}
								)
set_target_properties(pymass PROPERTIES COMPILE_FLAGS "-Dpymass_EXPORTS")
install(TARGETS pymass DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../python)
generate_export_header(pymass)


add_executable(testMZXML utils.cpp LCMS.cpp mzXMLParser.cpp main.cpp )
set_target_properties(testMZXML PROPERTIES COMPILE_FLAGS "-DPYMASS_STATIC_DEFINE")
target_include_directories(testMZXML PRIVATE ${TBB_INCLUDE_DIRS})
target_link_libraries(testMZXML ${EXPAT_LIBRARIES} base64
								debug     ${TBB_LIBRARIES_DEBUG}
								optimized ${TBB_LIBRARIES}
								)
