#include "CScoreProcessor/CNorm/CNorm.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdio>
#include <matio.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CNorm (L2) works for high snr inputs (q: 64)", "[norm][high][l2]") {

    constexpr unsigned q = 64;

    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data_w1_nofreq.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "data_input_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "data_input_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    mat_complex_split_t * in_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) in_data->Re, ((float *) in_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) in_data->Im, ((float *) in_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "norms_l2_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "norms_l2_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const vector<float> norms_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CNorm<q> * proc = new CNorm<q>();

    vector<float> results(norms_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        results[i] = proc->process(re_in[i], im_in[i]);
    }

    // Initialization Startup
    for (int64_t i = 0; i < int64_t(q - 1); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(norms_out[i] * norms_out[i] + 1, 1e-9f));
    }

    for (int64_t i = q - 1; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(std::sqrt(results[i]), Catch::Matchers::WithinRel(norms_out[i], 1e-9f));
        // printf("%ld\n", i);
    }

    delete proc;
}

TEST_CASE("CNorm (L2) works for low snr inputs (q: 64)", "[norm][low][l2]") {

    constexpr unsigned q = 64;

    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data_w1_nofreq.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "data_input_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "data_input_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    mat_complex_split_t * in_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) in_data->Re, ((float *) in_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) in_data->Im, ((float *) in_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "norms_l2_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "norms_l2_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const vector<float> norms_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CNorm<q> * proc = new CNorm<q>();

    vector<float> results(norms_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        results[i] = proc->process(re_in[i], im_in[i]);
    }

    // Initialization Startup
    for (int64_t i = 0; i < int64_t(q - 1); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(norms_out[i] * norms_out[i] + 1, 1e-5f));
    }

    for (int64_t i = q - 1; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(std::sqrt(results[i]), Catch::Matchers::WithinRel(norms_out[i], 5e-4f));
    }

    delete proc;
}
