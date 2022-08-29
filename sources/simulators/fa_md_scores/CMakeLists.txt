
find_package (OpenMP REQUIRED)
find_package (Boost COMPONENTS program_options REQUIRED)

add_executable (fa_md_simulation main.cpp)

target_link_directories (fa_md_simulation PRIVATE ${matio_INCLUDEDIR})
target_link_libraries (fa_md_simulation PRIVATE generator detectors Boost::program_options OpenMP::OpenMP_CXX ${matio} nbldpc_wrap)
