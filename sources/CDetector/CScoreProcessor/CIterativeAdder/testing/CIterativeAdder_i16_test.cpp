#include "CScoreProcessor/CIterativeAdder/CIterativeAdder.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <matio.h>
#include <matioCpp/matioCpp.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CIterativeAdder int16_t works for high snr inputs (q: 64, w:16, i:2)", "[iterativeadder][high][fixed]") {

    constexpr unsigned q     = 64;
    constexpr unsigned In_W  = 16;
    constexpr unsigned In_I  = 2;
    constexpr unsigned Out_W = In_W + 1;
    constexpr unsigned Out_I = In_I + 1;

    constexpr float in_scale_factor = float(1 << (In_W - In_I));
    constexpr float ot_scale_factor = 1.f / float(1 << (Out_W - Out_I));

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

    CIterativeAdder<q, int16_t> * proc = new CIterativeAdder<q, int16_t>();

    vector<float> re_results(re_out.size(), 0.f);
    vector<float> im_results(im_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(re_results.size()); i++) {
        const int16_t fx_re_in = int16_t(re_in[i] * in_scale_factor);
        const int16_t fx_im_in = int16_t(im_in[i] * in_scale_factor);

        int32_t fx_re_out, fx_im_out;

        proc->process(fx_re_in, fx_im_in, &fx_re_out, &fx_im_out);

        re_results[i] = float(fx_re_out) * ot_scale_factor;
        im_results[i] = float(fx_im_out) * ot_scale_factor;
    }

    for (int64_t i = 0; i < int64_t(re_results.size()); i++) {
        REQUIRE_THAT(re_results[i], Catch::Matchers::WithinRel(re_out[i], 1e-4f));
        REQUIRE_THAT(im_results[i], Catch::Matchers::WithinRel(im_out[i], 1e-4f));
    }

    delete proc;
}

TEST_CASE("CIterativeAdder int16_t works for low snr inputs (q: 64, w: 16, i:5)", "[iterativeadder][low][fixed]") {

    constexpr unsigned q     = 64;
    constexpr unsigned In_W  = 16;
    constexpr unsigned In_I  = 5;
    constexpr unsigned Out_W = In_W + 1;
    constexpr unsigned Out_I = In_I + 1;

    constexpr float in_scale_factor = float(1 << (In_W - In_I));
    constexpr float ot_scale_factor = 1.f / float(1 << (Out_W - Out_I));

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

    CIterativeAdder<q, int16_t> * proc = new CIterativeAdder<q, int16_t>();

    vector<float> re_results(re_out.size(), 0.f);
    vector<float> im_results(im_out.size(), 0.f);

    constexpr float margin = 12 * 1e-4; // -12 < Data input < 12

    for (int64_t i = 0; i < int64_t(re_results.size()); i++) {
        const int16_t fx_re_in = int16_t(re_in[i] * in_scale_factor);
        const int16_t fx_im_in = int16_t(im_in[i] * in_scale_factor);

        int32_t fx_re_out, fx_im_out;

        REQUIRE_THAT(float(fx_re_in), Catch::Matchers::WithinAbs(re_in[i] * in_scale_factor, margin * in_scale_factor));
        REQUIRE_THAT(float(fx_im_in), Catch::Matchers::WithinAbs(im_in[i] * in_scale_factor, margin * in_scale_factor));

        proc->process(fx_re_in, fx_im_in, &fx_re_out, &fx_im_out);

        re_results[i] = float(fx_re_out) * ot_scale_factor;
        im_results[i] = float(fx_im_out) * ot_scale_factor;
    }

    for (int64_t i = 0; i < int64_t(re_results.size()); i++) {
        REQUIRE_THAT(re_results[i], Catch::Matchers::WithinAbs(re_out[i], margin));
        REQUIRE_THAT(im_results[i], Catch::Matchers::WithinAbs(im_out[i], margin));
    }

    delete proc;
}
