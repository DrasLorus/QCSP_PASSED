add_library (nbldpc_wrap OBJECT init.cpp tool.cpp)
target_link_libraries (nbldpc_wrap PUBLIC ${LIB_MATH})
