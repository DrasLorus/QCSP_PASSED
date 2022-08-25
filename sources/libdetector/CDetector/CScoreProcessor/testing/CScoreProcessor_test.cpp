#include "../CScoreProcessor.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <matio.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CScoreProcessor works for high snr inputs (q: 64, N: 60)", "[scoreproc][high]") {

    constexpr unsigned q = 64;
    constexpr unsigned N = 60;

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

    tmp_mat = Mat_VarRead(data_file, "score_l2_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "score_l2_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const vector<float> score_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CScoreProcessor<N, q, float> * proc = new CScoreProcessor<N, q, float>(pn.data());

    vector<float> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        results[i] = proc->process(re_in[i], im_in[i]);
    }

    constexpr float margin = 480.f * 1e-4;

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(score_out[i], margin));
    }

    delete proc;
}

TEST_CASE("CScoreProcessor works for low snr inputs (q: 64, N: 60)", "[scoreproc][low]") {

    constexpr unsigned q = 64;
    constexpr unsigned N = 60;

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

    tmp_mat = Mat_VarRead(data_file, "score_l2_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "score_l2_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const vector<float> score_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CScoreProcessor<N, q, float> * proc = new CScoreProcessor<N, q, float>(pn.data());

    vector<float> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        results[i] = proc->process(re_in[i], im_in[i]);
    }

    constexpr float margin = 160.f * 5e-4f;

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(score_out[i], margin));
    }

    delete proc;
}
