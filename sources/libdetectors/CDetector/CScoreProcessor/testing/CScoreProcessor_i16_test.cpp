#include "../CScoreProcessor.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdio>
#include <matio.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CScoreProcessor int16_t works for high snr inputs (q: 64, N: 60)", "[scoreproc][high][fixed][.]") {

    constexpr unsigned q     = 64;
    constexpr unsigned N     = 60;
    constexpr unsigned In_W  = 16;
    constexpr unsigned In_I  = 7;
    constexpr unsigned Out_W = 24;
    constexpr unsigned Out_I = 23;

    constexpr float in_scale_factor = float(1U << (In_W - In_I));
    constexpr float ot_scale_factor = 1. / float(1 << (Out_W - Out_I));

    mat_t * parm_file = Mat_Open("../data/parameters_20210903.mat", MAT_ACC_RDONLY);
    if (not bool(parm_file)) {
        throw "../data/parameters_20210903.mat can't be opened.";
    }

    matvar_t * tmp_pn = Mat_VarRead(parm_file, "PN64");
    if (not bool(tmp_pn)) {
        throw "PN64 can't be loaded.";
    }

    const vector<float> pn((double *) tmp_pn->data, (double *) tmp_pn->data + tmp_pn->dims[1]);
    if (pn.size() != q) {
        throw "PN64 and q don't match.";
    }

    Mat_VarFree(tmp_pn);
    Mat_Close(parm_file);
    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data_w1_nofreq.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "data_input_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "data_input_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    mat_complex_split_t * cpx_ptr = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) cpx_ptr->Re, ((float *) cpx_ptr->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) cpx_ptr->Im, ((float *) cpx_ptr->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "score_sqr_l2_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "score_sqr_l2_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const vector<float> score_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CScoreProcessor<N, q, int16_t> * proc = new CScoreProcessor<N, q, int16_t>(pn.data());

    vector<float> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const int16_t fx_re_in = int16_t(re_in[i] * in_scale_factor);
        const int16_t fx_im_in = int16_t(im_in[i] * in_scale_factor);

        const uint64_t fx_score = proc->process_sqr(fx_re_in, fx_im_in);

        results[i] = float(fx_score) * ot_scale_factor;
    }

    for (int64_t i = q * N * 3 - 1; i < int64_t(results.size()); i += q * N * 5) { // Only maxima matters
        INFO("input no " << i);
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(score_out[i], 1.e-8f));
    }

    delete proc;
}

TEST_CASE("CScoreProcessor int16_t works for low snr inputs (q: 64, N: 60)", "[scoreproc][low][fixed]") {

    constexpr unsigned q     = 64;
    constexpr unsigned N     = 60;
    constexpr unsigned In_W  = 16;
    constexpr unsigned In_I  = 4;
    constexpr unsigned Out_W = 24;
    constexpr unsigned Out_I = 23;

    constexpr float in_scale_factor = float(1U << (In_W - In_I));
    constexpr float ot_scale_factor = 1.f / float(1 << (Out_W - Out_I));

    mat_t * parm_file = Mat_Open("../data/parameters_20210903.mat", MAT_ACC_RDONLY);
    if (not bool(parm_file)) {
        throw "../data/parameters_20210903.mat can't be opened.";
    }

    matvar_t * tmp_pn = Mat_VarRead(parm_file, "PN64");
    if (not bool(tmp_pn)) {
        throw "PN64 can't be loaded.";
    }

    const vector<float> pn((double *) tmp_pn->data, (double *) tmp_pn->data + tmp_pn->dims[1]);
    if (pn.size() != q) {
        throw "PN64 and q don't match.";
    }

    Mat_VarFree(tmp_pn);
    Mat_Close(parm_file);

    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data_w1_nofreq.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "data_input_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "data_input_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    mat_complex_split_t * cpx_ptr = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) cpx_ptr->Re, ((float *) cpx_ptr->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) cpx_ptr->Im, ((float *) cpx_ptr->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "score_sqr_l2_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "score_sqr_l2_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const vector<float> score_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CScoreProcessor<N, q, int16_t> * proc = new CScoreProcessor<N, q, int16_t>(pn.data());

    vector<float> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const int16_t fx_re_in = int16_t(re_in[i] * in_scale_factor);
        const int16_t fx_im_in = int16_t(im_in[i] * in_scale_factor);

        const uint64_t fx_score = proc->process_sqr(fx_re_in, fx_im_in);

        results[i] = float(fx_score) * ot_scale_factor;
    }

    for (int64_t i = int64_t(q - 1); i < int64_t(results.size()); i++) {
        if (score_out[i] > 330.f) { // Only focus on maxima
            INFO("input no " << i);
            REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(score_out[i], 50));
        }
    }

    delete proc;
}
