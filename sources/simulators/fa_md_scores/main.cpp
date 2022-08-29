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

constexpr int defines_line                    = __LINE__;
constexpr int p_omega_count                   = 9;
constexpr int p_omega_accepted[p_omega_count] = {1, 2, 3, 4, 5, 6, 7, 8, 9};
constexpr int q_count                         = 4;
constexpr int q_accepted[q_count]             = {64, 128, 256, 512};
constexpr int N_count                         = 1;
constexpr int N_accepted[1]                   = {60};

#define switch_create_p_omega(ptr, CNidx, CQidx, Oidx)                                                      \
    switch (Oidx) {                                                                                         \
        case 0:                                                                                             \
            ptr = new vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[0]>>(); \
            break;                                                                                          \
        case 1:                                                                                             \
            ptr = new vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[1]>>(); \
            break;                                                                                          \
        case 2:                                                                                             \
            ptr = new vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[2]>>(); \
            break;                                                                                          \
        case 3:                                                                                             \
            ptr = new vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[3]>>(); \
            break;                                                                                          \
        case 4:                                                                                             \
            ptr = new vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[4]>>(); \
            break;                                                                                          \
        case 5:                                                                                             \
            ptr = new vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[5]>>(); \
            break;                                                                                          \
        case 6:                                                                                             \
            ptr = new vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[6]>>(); \
            break;                                                                                          \
        case 7:                                                                                             \
            ptr = new vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[7]>>(); \
            break;                                                                                          \
        case 8:                                                                                             \
            ptr = new vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[8]>>(); \
            break;                                                                                          \
    }

#define switch_create_q(ptr, CNidx, Qidx, Oidx)         \
    switch (Qidx) {                                     \
        case 0:                                         \
            switch_create_p_omega(ptr, CNidx, 0, Oidx); \
            break;                                      \
        case 1:                                         \
            switch_create_p_omega(ptr, CNidx, 1, Oidx); \
            break;                                      \
        case 2:                                         \
            switch_create_p_omega(ptr, CNidx, 2, Oidx); \
            break;                                      \
        case 3:                                         \
            switch_create_p_omega(ptr, CNidx, 3, Oidx); \
            break;                                      \
    }

#define allocate_detector_vector(ptr, Nidx, Qidx, Oidx) \
    switch (Nidx) {                                     \
        case 0:                                         \
            switch_create_q(ptr, 0, Qidx, Oidx);        \
            break;                                      \
    }

#define switch_delete_p_omega(ptr, CNidx, CQidx, Oidx)                                                         \
    switch (Oidx) {                                                                                            \
        case 0:                                                                                                \
            delete (vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[0]>> *) ptr; \
            break;                                                                                             \
        case 1:                                                                                                \
            delete (vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[1]>> *) ptr; \
            break;                                                                                             \
        case 2:                                                                                                \
            delete (vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[2]>> *) ptr; \
            break;                                                                                             \
        case 3:                                                                                                \
            delete (vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[3]>> *) ptr; \
            break;                                                                                             \
        case 4:                                                                                                \
            delete (vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[4]>> *) ptr; \
            break;                                                                                             \
        case 5:                                                                                                \
            delete (vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[5]>> *) ptr; \
            break;                                                                                             \
        case 6:                                                                                                \
            delete (vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[6]>> *) ptr; \
            break;                                                                                             \
        case 7:                                                                                                \
            delete (vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[7]>> *) ptr; \
            break;                                                                                             \
        case 8:                                                                                                \
            delete (vector<CDetectorSerial<N_accepted[CNidx], q_accepted[CQidx], p_omega_accepted[8]>> *) ptr; \
            break;                                                                                             \
    }

#define switch_delete_q(ptr, CNidx, Qidx, Oidx)         \
    switch (Qidx) {                                     \
        case 0:                                         \
            switch_delete_p_omega(ptr, CNidx, 0, Oidx); \
            break;                                      \
        case 1:                                         \
            switch_delete_p_omega(ptr, CNidx, 1, Oidx); \
            break;                                      \
        case 2:                                         \
            switch_delete_p_omega(ptr, CNidx, 2, Oidx); \
            break;                                      \
        case 3:                                         \
            switch_delete_p_omega(ptr, CNidx, 3, Oidx); \
            break;                                      \
    }

#define delete_detector_vector(ptr, Nidx, Qidx, Oidx) \
    switch (Nidx) {                                   \
        case 0:                                       \
            switch_delete_q(ptr, 0, Qidx, Oidx);      \
            break;                                    \
    }

int main(int argc, char * argv[]) {

    po::options_description desc("Options");
    desc.add_options()(
        "help,h", "produce help message")(
        "snr,s", po::value<float>()->default_value(-10.f), "targeted signal-to-noise ratio")(
        "runs,n", po::value<int>()->default_value(10), "number of monte-carlo runs")(
        "p_delta,d", po::value<int>()->default_value(1), "number of score calculated every q chips (No effect with time-sliding)")(
        "p_omega,w", po::value<int>()->default_value(1), "number of rotation hypotheses")(
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

    const int p_omega = vm["p_omega"].as<int>();
    if (std::find(p_omega_accepted, p_omega_accepted + p_omega_count, p_omega) == (p_omega_accepted + p_omega_count)) {
        std::cerr << "\033[31m[Error]\033[0m p_omega value " << p_omega << " is not currently supported.\n"
                  << "You can try amending " << __FILE__ << " at line " << defines_line + 1 << "." << std::endl;
        return 1;
    }

    // const int p_delta = vm["p_delta"].as<int>();

    constexpr float rotation_span = two_pi_f;

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

    if (std::find(q_accepted, q_accepted + q_count, q) == (q_accepted + q_count)) {
        std::cerr << "\033[31m[Error]\033[0m q value " << q << " is not currently supported.\n"
                  << "Amend " << __FILE__ << " at line " << defines_line + 3 << "." << std::endl;
        return 1;
    }

    if (std::find(N_accepted, N_accepted + N_count, N) == (N_accepted + N_count)) {
        std::cerr << "\033[31m[Error]\033[0m N value " << N << " is not currently supported.\n"
                  << "Amend " << __FILE__ << " at line " << defines_line + 5 << "." << std::endl;
        return 1;
    }

    const string pn_path = vm["pn-file"].as<string>();
    mat_t *      pn_file = Mat_Open(pn_path.c_str(), MAT_ACC_RDONLY);
    if (!bool(pn_file)) {
        cerr << "\033[31m[Error]\033[0m Can't open pn mat file " << pn_path << "\n  "
             << strerror(errno) << endl;
        exit(EXIT_FAILURE);
    }

    const string pn_name = "PN" + std::to_string(q);
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

    const vector<float> pn(pn_data, pn_data + q);

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

    vector<CGenerator> generators;
    for (int t = 0; t < threads; t++) {
        generators.emplace_back(alist_path, pn, ovmod, snr, rotation_span);
    }

    void * detectors;

    allocate_detector_vector(detectors, 0, 0, 0);

    

    delete_detector_vector(detectors, 0, 0, 0);

    return EXIT_SUCCESS;
}
