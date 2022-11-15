find_package (Threads REQUIRED)
find_package (Boost COMPONENTS program_options REQUIRED)

configure_file (config.h.in ${CMAKE_CURRENT_SOURCE_DIR}/config.h)

add_executable (detector_benchmark_q64_w1 main.cpp)

target_include_directories (detector_benchmark_q64_w1 BEFORE PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries (detector_benchmark_q64_w1 PRIVATE detectors ${matio} Boost::program_options Threads::Threads)
