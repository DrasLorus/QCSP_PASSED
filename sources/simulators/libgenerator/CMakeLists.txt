add_subdirectory (CGenerator/nbldpc_functions)

add_library (generator STATIC CGenerator/CGenerator.cpp)
target_include_directories (generator PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (generator PUBLIC nbldpc_wrap)
