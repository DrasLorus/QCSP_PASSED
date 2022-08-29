#include "CDetector/CDetector.hpp"
#include "CGenerator/CGenerator.hpp"

#include <iostream>
#include <matio.h>

#include <boost/program_options.hpp>
#include <omp.h>
#include <string>

using std::string;
using std::vector;

namespace po = boost::program_options;

using namespace QCSP::StandaloneDetector;
using namespace QCSP::StandaloneDetector::Simulators;

using std::cerr;
using std::cout;
using std::endl;

enum sim_type_t {
    FA_MD,
    FA,
    MD
};

static constexpr int cst_q       = 64;
static constexpr int cst_N       = 60;
static constexpr int cst_p_omega = 4;

int main(int argc, char * argv[]) {
    static_assert(cst_p_omega > 0, "p_omega cannot be null.");
    static_assert(is_pow2(cst_q), "q must be a power of 2.");
    static_assert(cst_N == 60, "cst_N must be 60.");

    po::options_description desc("Options");
    desc.add_options()(
        "help,h", "produce help message")(
        "snr,s", po::value<float>()->default_value(-10.f), "targeted signal-to-noise ratio")(
        "runs,n", po::value<int>()->default_value(10), "number of monte-carlo runs")(
        "p_delta,d", po::value<int>()->default_value(1), "number of score calculated every q chips (No effect with time-sliding)")(
        "step-denominator", po::value<int>()->default_value(1), "step denominator, maximum error = pi / (2 * step)")(
        "pn-file", po::value<string>()->required(), "pn sequence mat file to use")(
        "ovmod-file", po::value<string>()->required(), "overmodulation sequence file to use")(
        "alist-file", po::value<string>()->required(), "alist file to use")(
        "threads,t", po::value<int>()->default_value(1), "number of threads to use")(
        "full-score", "log the full score instead of only the maximum one (WARNING: LARGE FILES AND HEAVY PERFORMANCE IMPACT)")(
        "output-file,o", po::value<string>()->default_value("score.mat"), "score file to write");

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
    const float snr = vm["snr"].as<float>();
    // const int   runs = vm["runs"].as<int>();

    const int threads = vm["threads"].as<int>();

    // const int p_delta = vm["p_delta"].as<int>();

    const string alist_path = vm["alist-file"].as<string>();

    FILE * alist_file = fopen(alist_path.c_str(), "r");
    if (!bool(alist_file)) {
        cerr << "\033[31m[Error]\033[0m Can't open alist file " << alist_path << "\n  "
             << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }

    int       N, M, q;
    const int nmq_read = fscanf(alist_file, "%d %d %d", &N, &M, &q);
    if (nmq_read != 3) {
        cerr << "\033[31m[Error]\033[0m Alist file " << alist_path << " seems ill-formed:\n"
             << "  Failed to read \"N M q\" from it." << endl;
        fclose(alist_file);
        exit(EXIT_FAILURE);
    }

    fclose(alist_file);

    if (q != cst_q) {
        std::cerr << "\033[31m[Error]\033[0m q value (" << q
                  << ") is not currently the one the program have been compiled for (" << cst_q << ").\n"
                  << "You may have to reconfigure and rebuild the project to fix that." << std::endl;
        return 1;
    }

    if (N != cst_N) {
        std::cerr << "\033[31m[Error]\033[0m N value (" << N
                  << ") is not currently the one the program have been compiled for (" << cst_N << ").\n"
                  << "You may have to reconfigure and rebuild the project to fix that." << std::endl;
        return 1;
    }

    const string pn_path = vm["pn-file"].as<string>();
    mat_t *      pn_file = Mat_Open(pn_path.c_str(), MAT_ACC_RDONLY);
    if (!bool(pn_file)) {
        cerr << "\033[31m[Error]\033[0m Can't open pn mat file " << pn_path << "\n  "
             << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }

    const string pn_name = "PN" + std::to_string(cst_q);
    matvar_t *   pn_var  = Mat_VarRead(pn_file, pn_name.c_str());
    if (!bool(pn_var)) {
        cerr << "\033[31m[Error]\033[0m Can't find pn mat variable " << pn_name << " in " << pn_path << "\n  "
             << strerror(errno) << endl;
        Mat_Close(pn_file);
        exit(EXIT_FAILURE);
    }

    const float * pn_data = (float *) pn_var->data;
    const size_t  pn_size = pn_var->dims[1];

    if (pn_size != size_t(q)) {
        cerr << "\033[31m[Error]\033[0m alist q (" << q << ") and size of "
             << pn_name << " in " << pn_path << " (" << pn_size << ") don't match." << endl;
        Mat_VarFree(pn_var);
        Mat_Close(pn_file);
        exit(EXIT_FAILURE);
    }

    const vector<float> pn(pn_data, pn_data + cst_q);

    Mat_VarFree(pn_var);
    Mat_Close(pn_file);

    const string ovmod_path = vm["ovmod-file"].as<string>();
    mat_t *      ovmod_file = Mat_Open(ovmod_path.c_str(), MAT_ACC_RDONLY);
    if (!bool(ovmod_file)) {
        cerr << "\033[31m[Error]\033[0m Can't open ovmod mat file " << ovmod_path << "\n  "
             << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }

    const string ovmod_name = "best_N";
    matvar_t *   ovmod_var  = Mat_VarRead(ovmod_file, ovmod_name.c_str());
    if (!bool(ovmod_var)) {
        cerr << "\033[31m[Error]\033[0m Can't find ovmod mat variable " << ovmod_name << " in " << ovmod_path << "\n  "
             << strerror(errno) << endl;
        Mat_Close(ovmod_file);
        exit(EXIT_FAILURE);
    }

    const double * ovmod_data = (double *) ovmod_var->data;
    const size_t   ovmod_size = ovmod_var->dims[1];

    if (ovmod_size != size_t(N)) {
        cerr << "\033[31m[Error]\033[0m alist N (" << N << ") and size of "
             << ovmod_name << " in " << ovmod_path << " (" << ovmod_size << ") don't match." << endl;
        Mat_VarFree(ovmod_var);
        Mat_Close(ovmod_file);
        exit(EXIT_FAILURE);
    }

    const vector<float> ovmod(ovmod_data, ovmod_data + N);

    Mat_VarFree(ovmod_var);
    Mat_Close(pn_file);

    const int step_denominator = vm["step-denominator"].as<int>();

    const float rotation_span = (cst_p_omega - 1) * two_pi_f / step_denominator;

    vector<CGenerator> generators;
    for (int t = 0; t < threads; t++) {
        generators.emplace_back(alist_path, pn, ovmod, snr, rotation_span);
    }

    vector<CDetectorSerial<cst_N, cst_q, cst_p_omega, float>> detectors;
    for (int t = 0; t < threads; t++) {
        detectors.emplace_back(pn.data(), 450, step_denominator);
    }
    cout << "Symbol rotation: " << (detectors[0].symbol_rotation() / pi_f) << " pi" << endl;
    cout << "Rotation span:    " << (rotation_span / pi_f) << " pi" << endl;

    return EXIT_SUCCESS;
}
