#ifndef _QCSP_PASSED_BENCHMARK_CONFIG_H_
#define _QCSP_PASSED_BENCHMARK_CONFIG_H_

#include <iostream>

/* #undef USE_UNISTD_API */
/* #undef USE_WINDOWS_API */

#define error_stream   std::cerr << "\033[31;1m[Error]\033[0m "
#define warning_stream std::cerr << "\033[33;1m[Warning]\033[0m "
#define info_stream    std::cerr << "\033[34;1m[INFO]\033[0m "

#ifndef NDEBUG
#define debug_print(info) cerr << "\033[35;1m[DEBUG]\033[0m " << info << endl
#else
#include <cassert>
#define debug_print(info) assert(true)
#endif

#define USE_MULT 0
#if defined(USE_MULT) && (USE_MULT == 1)
#define TDetector    CDetectorSerialMult
#define STR_DETECTOR "tsmult"
#else
#define TDetector    CDetectorSerial
#define STR_DETECTOR "ts"
#endif

#endif // _QCSP_PASSED_BENCHMARK_CONFIG_H_
