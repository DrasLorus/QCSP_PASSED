
add_subdirectory (nbldpc_functions)

add_library (generator STATIC CGenerator.cpp)
target_link_libraries (generator PUBLIC nbldpc_wrap)
