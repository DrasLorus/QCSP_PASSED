#include "../CNorm.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdint>
#include <cstdio>
#include <matio.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CNorm (L2) int16_t works for high snr inputs (q: 64)", "[norm][high][l2][fixed]") {

    constexpr unsigned q     = 64;
    constexpr unsigned In_W  = 16;
    constexpr unsigned In_I  = 4;
    constexpr unsigned Out_W = 33;
    constexpr unsigned Out_I = 9;

    constexpr float in_scale_factor = float(1U << (In_W - In_I));
    constexpr float ot_scale_factor = 1. / float(1 << (Out_W - Out_I));

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

    tmp_mat = Mat_VarRead(data_file, "norms_l2_sqr_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "norms_l2_sqr_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const vector<float> norms_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CNorm<q, int16_t> * proc = new CNorm<q, int16_t>();

    vector<float> results(norms_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const int16_t fx_re_in = int16_t(re_in[i] * in_scale_factor);
        const int16_t fx_im_in = int16_t(im_in[i] * in_scale_factor);

        const uint64_t fx_norm = proc->process_sqr(fx_re_in, fx_im_in);

        results[i] = float(fx_norm) * ot_scale_factor;
    }

    // Initialization Startup
    for (int64_t i = 0; i < int64_t(q); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(norms_out[i] + ot_scale_factor, 1e-7f));
    }

    for (int64_t i = q; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(norms_out[i], 1e-9f));
        // printf("%ld\n", i);
    }

    delete proc;
}

TEST_CASE("CNorm (L2) int16_t works for low snr inputs (q: 64)", "[norm][low][l2][fixed]") {

    constexpr unsigned q     = 64;
    constexpr unsigned In_W  = 16;
    constexpr unsigned In_I  = 6;
    constexpr unsigned Out_W = 33;
    constexpr unsigned Out_I = 13;

    constexpr float in_scale_factor = float(1U << (In_W - In_I));
    constexpr float ot_scale_factor = 1. / float(1 << (Out_W - Out_I));

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

    tmp_mat = Mat_VarRead(data_file, "norms_l2_sqr_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "norms_l2_sqr_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const vector<float> norms_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CNorm<q, int16_t> * proc = new CNorm<q, int16_t>();

    vector<float> results(norms_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const int16_t fx_re_in = int16_t(re_in[i] * in_scale_factor);
        const int16_t fx_im_in = int16_t(im_in[i] * in_scale_factor);

        const uint64_t fx_norm = proc->process_sqr(fx_re_in, fx_im_in);

        results[i] = float(fx_norm) * ot_scale_factor;
    }

    // Initialization Startup
    for (int64_t i = 0; i < int64_t(q); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(norms_out[i] + ot_scale_factor, .5f));
    }

    for (int64_t i = q; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(norms_out[i], 1e-3f));
        // printf("%ld\n", i);
    }

    delete proc;
}
