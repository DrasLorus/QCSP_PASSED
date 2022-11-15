#include <boost/program_options.hpp>
#include <matio.h>
#include <string>

#include "./config.h"

#include "CDetector/CDetector.hpp"
#include "Miscellanous/misc.hpp"

namespace po = boost::program_options;

using std::string;
using std::vector;

constexpr unsigned N       = 60;
constexpr unsigned q       = 64;
constexpr unsigned p_omega = 5;

typedef QCSP::StandaloneDetector::CDetectorSerial<N, q, p_omega, float, true> detector_t;
using state_t = detector_t::state_t;

int main(int argc, char * argv[]) {

    po::options_description desc("Options");
    desc.add_options()(
        "help,h", "produce help message")(
        "pn-file", po::value<string>()->required(), "pn sequence mat file to use")(
        "threshold,H", po::value<float>()->default_value(450.f), "threshold to use in detection (no effect if complete is not specified)")(
        "complete", "simulate a complete detector instead of a partially inhibited one (WARNING: THRESHOLD MATTERS)");

    po::variables_map vm;
    try {
        po::store(po::command_line_parser(argc, argv)
                      .options(desc)
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

    for (int i = 0; i < int(1e8); i++) {
        detector->process(0.56f, -0.56f, &state);
    }
}
