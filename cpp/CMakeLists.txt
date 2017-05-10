include(GenerateExportHeader)
add_library(pymass SHARED utils.cpp base64.cpp LCMS.cpp mzXMLParser.cpp)
target_link_libraries(pymass ${EXPAT_LIBRARIES} base64)
set_target_properties(pymass PROPERTIES COMPILE_FLAGS "-Dpymass_EXPORTS")
install(TARGETS pymass DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/../python)
generate_export_header(pymass)


add_executable(testMZXML utils.cpp base64.cpp LCMS.cpp mzXMLParser.cpp main.cpp )
set_target_properties(testMZXML PROPERTIES COMPILE_FLAGS "-DPYMASS_STATIC_DEFINE")
target_link_libraries(testMZXML ${EXPAT_LIBRARIES} base64 gsort)
