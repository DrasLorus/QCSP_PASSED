#include "../CCorrAbsMax.hpp"
#include "Miscellanous/misc.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <cstdio>
#include <matio.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CCorrAbsMax int16_t works for high snr inputs (q: 64)", "[corrabsmax][high][fixed]") {
    constexpr size_t q     = 64;
    constexpr size_t p     = pow2_log2<q>();
    constexpr size_t In_W  = 17;
    constexpr size_t In_I  = 8; // High SNR => High signal power !
    constexpr size_t Out_W = 2 * (In_W + p + 1) + 1 - (p + 1) - (In_W + 1);
    constexpr size_t Out_I = 2 * (In_I + p + 1) + 1 - (p + 1) - 1;

    constexpr float  in_scale_factor = constexpr_pow<float, int64_t(In_W - In_I)>(2.);
    constexpr double ot_scale_factor = 1. / constexpr_pow<double, int64_t(Out_W - Out_I)>(2.);

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

    matvar_t * tmp_mat = Mat_VarRead(data_file, "iter_fcts_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "iter_fcts_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    mat_complex_split_t * in_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) in_data->Re, ((float *) in_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) in_data->Im, ((float *) in_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "cabs_max_sqr_raw_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "cabs_max_sqr_raw_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const vector<float> cabs_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CCorrAbsMax<q, int16_t> * proc = new CCorrAbsMax<q, int16_t>(pn.data());

    vector<float> results(cabs_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const int32_t fx_re_in = int32_t(re_in[i] * in_scale_factor);
        const int32_t fx_im_in = int32_t(im_in[i] * in_scale_factor);

        const uint64_t fx_out = proc->process(fx_re_in, fx_im_in);

        results[i] = float(double(fx_out) * ot_scale_factor);
    }

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        INFO("input no " << i);
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(cabs_out[i], 1e-9));
    }

    delete proc;
}

TEST_CASE("CCorrAbsMax int16_t works for low snr inputs (q: 64)", "[corrabsmax][low][fixed]") {

    constexpr unsigned q     = 64;
    constexpr uint     p     = pow2_log2<q>();
    constexpr unsigned In_W  = 17;
    constexpr unsigned In_I  = 5;
    constexpr uint     Out_W = 2 * (In_W + p + 1) + 1 - (p + 1) - (In_W + 1);
    constexpr uint     Out_I = 2 * (In_I + p + 1) + 1 - (p + 1) - 1;

    constexpr float  in_scale_factor = constexpr_pow<float, int64_t(In_W - In_I)>(2.);
    constexpr double ot_scale_factor = 1. / constexpr_pow<double, int64_t(Out_W - Out_I)>(2.);

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

    matvar_t * tmp_mat = Mat_VarRead(data_file, "iter_fcts_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "iter_fcts_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    mat_complex_split_t * in_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) in_data->Re, ((float *) in_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) in_data->Im, ((float *) in_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "cabs_max_sqr_raw_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "cabs_max_sqr_raw_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const vector<float> cabs_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CCorrAbsMax<q, int16_t> * proc = new CCorrAbsMax<q, int16_t>(pn.data());

    vector<float> results(cabs_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(64 * 60 * 5); i++) {
        const int32_t fx_re_in = int32_t(re_in[i] * in_scale_factor);
        const int32_t fx_im_in = int32_t(im_in[i] * in_scale_factor);

        const uint64_t fx_out = proc->process(fx_re_in, fx_im_in);

        results[i] = float(double(fx_out) * ot_scale_factor);
    }

    for (int64_t i = 0; i < int64_t(64 * 60 * 5); i++) {
        INFO("input no " << i);
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(cabs_out[i], 3e-2f));
    }

    delete proc;
}
