option (ENABLE_FAMD_MULT "Wether or not you want to instantiate plain multipliers. No benefit." OFF)
set (USE_MULT ${ENABLE_FAMD_MULT} CACHE INTERNAL "" FORCE)

add_subdirectory (latency_throughput)
add_subdirectory (libgenerator)
add_subdirectory (fa_md_scores)
