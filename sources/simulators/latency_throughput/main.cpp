#include <atomic>
#include <boost/program_options.hpp>
#include <boost/program_options/positional_options.hpp>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <matio.h>
#include <ostream>
#include <string>

#include "./config.h"

#include "CDetector/CDetector.hpp"
#include "Miscellanous/misc.hpp"

namespace po = boost::program_options;

using std::atomic;
using std::string;
using std::vector;

constexpr unsigned N       = 60;
constexpr unsigned q       = 64;
constexpr unsigned p_omega = 5;

typedef QCSP::StandaloneDetector::CDetectorSerial<N, q, p_omega, float, true> detector_t;
using state_t = detector_t::state_t;

static atomic<bool> stop_signal_not_called {true};

#if !defined(USE_WINDOWS_API) && defined(USE_UNISTD_API)
#include <csignal>

void sig_int_handler(int) {
    std::cerr << "\n  Stop signal received.\n  Stopping..." << std::endl;
    stop_signal_not_called = false;
}
#else
#include <windows.h>

BOOL WINAPI CtrlHandler(DWORD fdwCtrlType) {
    switch (fdwCtrlType) {
        // Handle the CTRL-C signal.
        case CTRL_C_EVENT:
            std::cerr << "\n\n * Stop signal received.\n * Stopping...\n"
                      << std::endl;
            stop_signal_not_called = false;
            return TRUE;
        default:
            return FALSE;
    }
}
#endif

int main(int argc, char * argv[]) {

    po::options_description desc("Options");
    desc.add_options()(
        "help,h", "produce help message")(
        "pn-file,p", po::value<string>()->required(), "pn sequence mat file to use")(
        "number-detections,n", po::value<int>()->default_value(int(-1)), "number of detections to perform.")(
        "threshold,H", po::value<float>()->default_value(450.f), "threshold to use in detection.");

    po::positional_options_description pod;
    pod.add("number-detections", 1);

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(desc)
                      .positional(pod)
                      .run(),
                  vm);

        if (vm.count("help")) {
            std::cout << "Usage: " << argv[0] << " [options]\n"
                      << desc << std::endl;
            return 1;
        }

        po::notify(vm);
    } catch (std::exception & e) {
        std::cerr << e.what() << std::endl;
        return 1;
    }

    const string pn_path = vm["pn-file"].as<string>();
    mat_t *      pn_file = Mat_Open(pn_path.c_str(), MAT_ACC_RDONLY);
    if (!bool(pn_file)) {
        error_stream << "Can't open pn mat file " << pn_path << std::endl;
        exit(EXIT_FAILURE);
    }

    const string pn_name = "PN" + std::to_string(q);
    matvar_t *   pn_var  = Mat_VarRead(pn_file, pn_name.c_str());
    if (!bool(pn_var)) {
        error_stream << "Can't find pn mat variable " << pn_name << " in " << pn_path << std::endl;
        Mat_Close(pn_file);
        exit(EXIT_FAILURE);
    }

    const float * pn_data = (float *) pn_var->data;
    const size_t  pn_size = pn_var->dims[1];

    if (pn_size != size_t(q)) {
        error_stream << "alist q (" << q << ") and size of "
                     << pn_name << " in " << pn_path << " (" << pn_size << ") don't match." << std::endl;
        Mat_VarFree(pn_var);
        Mat_Close(pn_file);
        exit(EXIT_FAILURE);
    }

    const vector<float> pn(pn_data, pn_data + q);

    Mat_VarFree(pn_var);
    Mat_Close(pn_file);

    auto conditional_n = [](int x) {
        return x >= 0 ? x : INT_MAX;
    };

    const int num_detections = conditional_n(vm.at("number-detections").as<int>());

#if !defined(USE_WINDOWS_API) && defined(USE_UNISTD_API)
    struct sigaction act;
    sigemptyset(&act.sa_mask);
    act.sa_handler = &sig_int_handler;
    act.sa_flags   = 0;

    if (sigaction(SIGINT, &act, nullptr)) {
        perror("'sigaction' failed.");
        exit(EXIT_FAILURE);
    }
#else
    if (not SetConsoleCtrlHandler(CtrlHandler, TRUE)) {
        std::cerr << "\nERROR: Could not set control handler" << std::endl;
        exit(EXIT_FAILURE);
    }
#endif

    detector_t * detector = new detector_t(pn.data(), 450.f, 0);

    state_t state {};
    state.frame_detected      = false;
    state.max_found           = false;
    state.max_score           = 0.f;
    state.frequency_offset    = 0.f;
    state.frequency_index     = 0;
    state.chip_from_max       = 0;
    state.chip_since_last_det = 0;
    std::memset(state.scores, 0, p_omega * sizeof(float));

    int i = 0;
    while ((i++ <= num_detections) && stop_signal_not_called) {
        detector->process(0.56f, -0.56f, &state);
    }
}
