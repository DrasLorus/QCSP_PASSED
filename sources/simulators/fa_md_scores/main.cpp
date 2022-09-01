#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>

#include <boost/program_options.hpp>
#include <libgen.h>
#include <string>
#include <unistd.h>

#include "./config.h"

#if defined(HAVE_UNISTD_H) && (HAVE_UNISTD_H == 1)
#include <libgen.h>
#include <sys/wait.h>
#include <unistd.h>
#endif

using std::string;

namespace po = boost::program_options;

using std::cerr;
// using std::cout;
using std::endl;

#define error_stream   cerr << "\033[31;1m[Error]\033[0m "
#define warning_stream cerr << "\033[33;1m[Warning]\033[0m "
#define info_stream    cerr << "\033[34;1m[INFO]\033[0m "

#ifndef NDEBUG
#define debug_print(info) cerr << "\033[35;1m[DEBUG]\033[0m " << info
#else
#define debug_print(info) assert(true)
#endif

int main(int argc, char * argv[]) {

    po::options_description desc("Options");
    desc.add_options()(
        "help,h", "produce help message")(
        "snr,s", po::value<float>()->default_value(-10.f), "targeted signal-to-noise ratio")(
        "runs,n", po::value<int>()->default_value(10), "number of monte-carlo runs")(
        "pn-file", po::value<string>()->required(), "pn sequence mat file to use")(
        "ovmod-file", po::value<string>()->required(), "overmodulation sequence file to use")(
        "alist-file", po::value<string>()->required(), "alist file to use")(
        "threads,t", po::value<int>()->default_value(1), "number of threads to use")(
        "delta,d", po::value<int>()->default_value(1), "value of p_delta to use (ignored for time sliding)")(
        "omega,w", po::value<int>()->required(), "value of p_omega to use")(
        "raw", "disable normalization")(
        "mult", "use mult-based detector")(
        "full-score", "log the full score instead of only the maximum one (WARNING: LARGE FILES AND HEAVY PERFORMANCE IMPACT)")(
        "output-file,o", po::value<string>()->default_value(DEFAULT_SCORE_FILE), "score file to write");

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

    const float  snr        = vm["snr"].as<float>();
    const int    runs       = vm["runs"].as<int>();
    const string pn_path    = vm["pn-file"].as<string>();
    const string ovmod_path = vm["ovmod-file"].as<string>();
    const string alist_path = vm["alist-file"].as<string>();
    const int    threads    = vm["threads"].as<int>();
    // const int    p_delta     = vm["delta"].as<int>();
    const int    p_omega     = vm["omega"].as<int>();
    const bool   normed      = !vm.count("raw");
    const bool   full_score  = vm.count("full-score");
    const string output_path = vm["output-file"].as<string>();

    FILE * alist = fopen(alist_path.c_str(), "r");
    if (!bool(alist)) {
        error_stream << "Failed to open alist " << alist_path << "." << endl;
        exit(EXIT_FAILURE);
    }

    int       N, q;
    const int result_scan = fscanf(alist, "%d %*d %d", &N, &q);
    if (result_scan != 2) {
        error_stream << "Failed to read N and q correctly in " << alist_path << ". Please, check its format." << endl;
        fclose(alist);
        exit(EXIT_FAILURE);
    }
    fclose(alist);

    const bool use_mult = vm.count("mult");

    std::ostringstream os;
    os << "fa_md_simulation_N" << N
       << "_q" << q
       << "_w" << p_omega
       << "_" << int(normed)
       << (use_mult ? "_mult" : "");

    const string program_name = os.str();

#if defined(_POSIX_VERSION)
    char * abs_path    = realpath(argv[0], nullptr);
    char * program_dir = dirname(abs_path);

    const string program_path = string(program_dir) + "/fa_md_simulation/" + program_name;

    free(abs_path);

    const string snr_str  = std::to_string(snr);
    const string runs_str = std::to_string(runs);

    char *       abs_pn = realpath(pn_path.c_str(), nullptr);
    const string pn_abs_path((bool(abs_pn) ? abs_pn : "pn-file-NOT-FOUND"));
    free(abs_pn);

    char *       abs_ovmod = realpath(ovmod_path.c_str(), nullptr);
    const string ovmod_abs_path((bool(abs_ovmod) ? abs_ovmod : "ovmod-file-NOT-FOUND"));
    free(abs_ovmod);

    char *       abs_alist = realpath(alist_path.c_str(), nullptr);
    const string alist_abs_path(abs_alist);
    free(abs_alist);

    const string threads_str = std::to_string(threads);

    const pid_t pid = fork();
    if (pid == 0) {
        const int result = execl(
            program_path.c_str(),
            program_name.c_str(),
            "--snr", snr_str.c_str(),
            "--runs", runs_str.c_str(),
            "--pn-file", pn_abs_path.c_str(),
            "--ovmod-file", ovmod_abs_path.c_str(),
            "--alist-file", alist_abs_path.c_str(),
            "--threads", threads_str.c_str(),
            (full_score ? "--full-score" : ""),
            "--output-file", output_path.c_str(),
            (char *) NULL);

        if (result != 0) {
            error_stream << "Failed to execute " << program_name << ": "
                         << (errno == ENOENT ? "Unsupported N, q or p_omega. Check CMake configuration." : strerror(errno))
                         << endl;
            exit(EXIT_FAILURE);
        }
    } else if (pid > 0) {
        int result;

        const pid_t child_pid = wait(&result);
        (void) child_pid; // Currently unused

        if (result != EXIT_SUCCESS) {
            error_stream << "Unknown error occured in " << program_name << " (" << (int) child_pid << ") with code " << result << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        perror("Failed to fork the process.");
        exit(EXIT_FAILURE);
    }
#endif
    return EXIT_SUCCESS;
}
