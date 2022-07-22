#include "CScoreProcessor/CIterativeAdder/CIterativeAdder.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <matio.h>
#include <matioCpp/matioCpp.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CIterativeAdder works for high snr inputs (q: 64)", "[iterativeadder][high]") {

    constexpr unsigned q = 64;

    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data_w1_nofreq.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "data_input_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "data_input_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const mat_complex_split_t * in_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) in_data->Re, ((float *) in_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) in_data->Im, ((float *) in_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "iter_fcts_infdB_w1_q64_N60_0_n10");
    if (not bool(tmp_mat)) {
        throw "iter_fcts_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const mat_complex_split_t * out_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_out((float *) out_data->Re, ((float *) out_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_out((float *) out_data->Im, ((float *) out_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CIterativeAdder<q> * proc = new CIterativeAdder<q>();

    vector<float> re_results(re_out.size(), 0.f);
    vector<float> im_results(im_out.size(), 0.f);

    float * p_re_out = re_results.data();
    float * p_im_out = im_results.data();
    for (int64_t i = 0; i < int64_t(re_results.size()); i++) {
        proc->process(re_in[i], im_in[i], p_re_out++, p_im_out++);
    }

    for (int64_t i = 0; i < int64_t(re_results.size()); i++) {
        REQUIRE_THAT(re_results[i], Catch::Matchers::WithinRel(re_out[i], 1e-9f));
        REQUIRE_THAT(im_results[i], Catch::Matchers::WithinRel(im_out[i], 1e-9f));
    }

    delete proc;
}

TEST_CASE("CIterativeAdder works for low snr inputs (q: 64)", "[iterativeadder][low]") {

    constexpr unsigned q = 64;

    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data_w1_nofreq.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "data_input_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "data_input_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const mat_complex_split_t * in_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) in_data->Re, ((float *) in_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) in_data->Im, ((float *) in_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    tmp_mat = Mat_VarRead(data_file, "iter_fcts_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "iter_fcts_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const mat_complex_split_t * out_data = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_out((float *) out_data->Re, ((float *) out_data->Re) + tmp_mat->dims[0]);
    const vector<float> im_out((float *) out_data->Im, ((float *) out_data->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    Mat_Close(data_file);

    CIterativeAdder<q> * proc = new CIterativeAdder<q>();

    vector<float> re_results(re_out.size(), 0.f);
    vector<float> im_results(im_out.size(), 0.f);

    float * p_re_out = re_results.data();
    float * p_im_out = im_results.data();
    for (int64_t i = 0; i < int64_t(re_results.size()); i++) {
        proc->process(re_in[i], im_in[i], p_re_out++, p_im_out++);
    }

    constexpr float margin = 12e-6; // -12 < Data input < 12

    for (int64_t i = 0; i < int64_t(re_results.size()); i++) {
        REQUIRE_THAT(re_results[i], Catch::Matchers::WithinAbs(re_out[i], margin));
        REQUIRE_THAT(im_results[i], Catch::Matchers::WithinAbs(im_out[i], margin));
    }

    delete proc;
}
