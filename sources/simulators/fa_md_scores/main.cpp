#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>

#include <boost/program_options.hpp>

#include "./config.h"

#if defined(HAVE_UNISTD_H) && (HAVE_UNISTD_H == 1)
#include <libgen.h>
#include <sys/wait.h>
#include <unistd.h>
#elif defined(_MSC_VER)
#include <windows.h>

#include <pathcch.h>
#include <processthreadsapi.h>
#endif

using std::string;

namespace po = boost::program_options;

using std::cerr;
// using std::cout;
using std::endl;

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
        "threshold,H", po::value<float>(), "threshold to use in detection (no effect if complete is not specified)")(
        "rotation-span", po::value<float>(), "rotation span to use. For a value X, rotation will be in the interval [-X/2, X/2]")(
        "step-numerator", po::value<unsigned>(), "step numerator to use")(
        "step-denominator", po::value<unsigned>(), "step denominator to use")(
        "complete", "simulate a complete detector instead of a partially inhibited one (WARNING: THRESHOLD MATTERS)")(
        "no-fa", "disable False Alarm scenario")(
        "no-md", "disable Miss Detection scenario")(
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

    const float    snr         = vm["snr"].as<float>();
    const int      runs        = vm["runs"].as<int>();
    const string   pn_path     = vm["pn-file"].as<string>();
    const string   ovmod_path  = vm["ovmod-file"].as<string>();
    const string   alist_path  = vm["alist-file"].as<string>();
    const int      threads     = vm["threads"].as<int>();
    const int      p_omega     = vm["omega"].as<int>();
    const bool     have_thres  = vm.count("threshold");
    const float    threshold   = have_thres ? vm["threshold"].as<float>() : 0.f;
    const bool     have_span   = vm.count("rotation-span");
    const float    rot_span    = have_span ? vm["rotation-span"].as<float>() : 0.f;
    const bool     have_num    = vm.count("step-numerator");
    const unsigned step_num    = have_num ? vm["step-numerator"].as<unsigned>() : 0U;
    const bool     have_step   = vm.count("step-denominator");
    const unsigned step_den    = have_step ? vm["step-denominator"].as<unsigned>() : 0U;
    const bool     complete    = vm.count("complete");
    const bool     no_fa       = vm.count("no-fa");
    const bool     no_md       = vm.count("no-md");
    const bool     normed      = !vm.count("raw");
    const bool     full_score  = vm.count("full-score");
    const string   output_path = vm["output-file"].as<string>();

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

    const string threads_str   = std::to_string(threads);
    const string threshold_str = std::to_string(threshold);
    const string span_str      = std::to_string(rot_span);
    const string num_str       = std::to_string(step_num);
    const string step_str      = std::to_string(step_den);

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
            "--output-file", output_path.c_str(),
            (have_thres ? "--threshold" : ""), (have_thres ? threshold_str.c_str() : ""),
            (have_span ? "--rotation-span" : ""), (have_span ? span_str.c_str() : ""),
            (have_num ? "--step-numerator" : ""), (have_num ? num_str.c_str() : ""),
            (have_step ? "--step-denominator" : ""), (have_step ? step_str.c_str() : ""),
            (complete ? "--complete" : ""),
            (no_fa ? "--no-fa" : ""),
            (no_md ? "--no-md" : ""),
            (full_score ? "--full-score" : ""),
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
            error_stream << "Unknown error occurred in " << program_name
                         << " (" << (int) child_pid << ") with code " << result << endl;
            exit(EXIT_FAILURE);
        }
    } else {
        perror("Failed to fork the process.");
        exit(EXIT_FAILURE);
    }
#elif defined(_MSC_VER)

    char *    abs_path    = _fullpath(nullptr, argv[0], 256);
    wchar_t * program_dir = new wchar_t[strlen(abs_path)];
    std::copy(abs_path, abs_path + strlen(abs_path), program_dir);

    const long int r = PathCchRemoveFileSpec(program_dir, wcslen(program_dir));
    if ((r != 0) && (r != 1)) {
        error_stream << "Path is ill-formed at line " << __LINE__ << " (code " << r << ")." << endl;
    }

    const string program_path = string(program_dir, program_dir + wcslen(program_dir)) + "\\fa_md_simulation\\" + program_name + ".exe";

    free(abs_path);
    delete[] program_dir;

    const string snr_str  = std::to_string(snr);
    const string runs_str = std::to_string(runs);

    char *       abs_pn = _fullpath(nullptr, pn_path.c_str(), pn_path.size());
    const string pn_abs_path((bool(abs_pn) ? abs_pn : "pn-file-NOT-FOUND"));
    free(abs_pn);

    char *       abs_ovmod = _fullpath(nullptr, ovmod_path.c_str(), ovmod_path.size());
    const string ovmod_abs_path((bool(abs_ovmod) ? abs_ovmod : "ovmod-file-NOT-FOUND"));
    free(abs_ovmod);

    char *       abs_alist = _fullpath(nullptr, alist_path.c_str(), alist_path.size());
    const string alist_abs_path(abs_alist);
    free(abs_alist);

    const string threads_str   = std::to_string(threads);
    const string threshold_str = std::to_string(threshold);
    const string span_str      = std::to_string(rot_span);
    const string num_str       = std::to_string(step_num);
    const string step_str      = std::to_string(step_den);

    std::ostringstream os_cmd;
    os_cmd << program_name << ".exe"
           << " --snr " << snr_str
           << " --runs " << runs_str
           << " --pn-file " << pn_abs_path
           << " --ovmod-file " << ovmod_abs_path
           << " --alist-file " << alist_abs_path
           << " --threads " << threads_str
           << " --output-file " << output_path
           << (have_thres ? " --threshold " : "") << (have_thres ? threshold_str : "")
           << (have_span ? " --rotation-span " : "") << (have_span ? span_str : "")
           << (have_num ? " --step-numerator " : "") << (have_num ? num_str : "")
           << (have_step ? " --step-denominator " : "") << (have_step ? step_str : "")
           << (complete ? " --complete " : "")
           << (no_fa ? " --no-fa " : "")
           << (no_md ? " --no-md " : "")
           << (full_score ? " --full-score " : "");
    const string cmdline_str = os_cmd.str();

    char * cmdline = new char[cmdline_str.size() + 1];
    for (int i = 0; i < cmdline_str.size(); i++) {
        cmdline[i] = cmdline_str[i];
    }
    cmdline[cmdline_str.size()] = (char) NULL;

    STARTUPINFOA        si;
    PROCESS_INFORMATION pi;

    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));

    const int success = CreateProcessA(
        program_path.c_str(),
        cmdline,
        nullptr,
        nullptr,
        TRUE,
        0,
        nullptr,
        nullptr,
        &si,
        &pi);

    if (!success) {
        error_stream << "CreateProcess failed for " << program_path << " (code " << GetLastError() << ")." << endl;
        delete[] cmdline;
        exit(EXIT_FAILURE);
    }

    // Wait until child process exits.
    WaitForSingleObject(pi.hProcess, INFINITE);

    // Close process and thread handles.
    CloseHandle(pi.hProcess);
    CloseHandle(pi.hThread);

    delete[] cmdline;
#endif
    return EXIT_SUCCESS;
}
