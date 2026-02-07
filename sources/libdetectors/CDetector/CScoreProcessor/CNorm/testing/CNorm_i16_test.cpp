#include "../CNorm.hpp"
#include "Miscellanous/metatypes.hpp"
#include "libdetectors/CDetector/CScoreProcessor/CNorm/CNormFP.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <iostream>
#include <matio.h>
#include <memory>
#include <numeric>
#include <vector>

using namespace QCSP::StandaloneDetector;
using namespace QCSP;

using std::vector;

static void error_compute(double   expected,
                          double   result,
                          double & total_sig,
                          double & total_qfx_noise,
                          double & total_rel_error,
                          double & total_sqr_error,
                          double & total_abs_error,
                          double & total_log_error) {
    const double qfx_noise = expected - result;
    const double rel_error = qfx_noise / std::abs(expected);
    total_sig += expected;
    total_qfx_noise += qfx_noise;
    total_rel_error += rel_error;
    total_sqr_error += std::pow(qfx_noise, 2);
    total_abs_error += std::abs(qfx_noise);
    total_log_error += std::pow(std::log10(expected + 1) - std::log10(result + 1), 2);
}

template <uint64_t In_W>
static typename saturation<32, 32 - In_W, metatypes::SIGNED>::out_t
input_stage(float value, float in_scale_factor) {
    const int32_t floored_v = (int32_t) std::floor(value * in_scale_factor);
    return saturation<32, 32 - In_W, metatypes::SIGNED>::apply(floored_v);
};

TEST_CASE("CNorm (L2) int16_t works for high snr inputs (q: 64)", "[norm][high][l2][fixed]") {

    constexpr unsigned q               = 64;
    constexpr int64_t  In_W            = 16;
    constexpr int64_t  In_I            = 6;
    constexpr int64_t  Out_W           = 16;
    constexpr int64_t  Out_I           = 12;
    constexpr float    in_scale_factor = constexpr_pow<float, In_W - In_I>(2);
    constexpr float    ot_scale_factor = constexpr_pow<float, Out_I - Out_W>(2);

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

    auto proc = std::make_unique<CNorm<q, intxx_t<In_W, metatypes::FROM_BITS>, In_W, In_I>>();

    using input_t = typename decltype(proc)::element_type::input_t;
    using out_t   = typename decltype(proc)::element_type::out_t;

    vector<float> results(norms_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const input_t fx_re_in = input_stage<In_W>(re_in[i], in_scale_factor);
        const input_t fx_im_in = input_stage<In_W>(im_in[i], in_scale_factor);

        const out_t fx_norm = proc->process_sqr(fx_re_in, fx_im_in);

        results[i] = float(fx_norm) * ot_scale_factor;
    }

    // Initialization Startup
    for (int64_t i = 0; i < int64_t(q - 1); i++) {
        constexpr float initial_out = float(decltype(proc)::element_type::initial_out * ot_scale_factor);
        INFO("input no " << i);
        REQUIRE_THAT(results[i] - initial_out, Catch::Matchers::WithinAbs(norms_out[i], ot_scale_factor));
    }

    for (int64_t i = q; i < int64_t(results.size()) - 1; i++) {
        INFO("input no " << i);
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(norms_out[i], ot_scale_factor));
    }
}

TEST_CASE("CNorm (L2) int16_t works for low snr inputs (q: 64)", "[norm][low][l2][fixed]") {

    constexpr unsigned q     = 64;
    constexpr int64_t  In_W  = 16;
    constexpr int64_t  In_I  = 6;
    constexpr int64_t  Out_W = 16;
    constexpr int64_t  Out_I = 12;

    constexpr float in_scale_factor = constexpr_pow<float, In_W - In_I>(2);
    // constexpr float ot_scale_factor = 1. / float(1 << (Out_W - Out_I));
    constexpr float ot_scale_factor = constexpr_pow<float, Out_I - Out_W>(2);

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

    auto proc = std::make_unique<CNorm<q, intxx_t<In_W, metatypes::FROM_BITS>, In_W, In_I>>();

    using input_t = typename decltype(proc)::element_type::input_t;
    using out_t   = typename decltype(proc)::element_type::out_t;

    vector<float> results(norms_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const input_t fx_re_in = input_stage<In_W>(re_in[i], in_scale_factor);
        const input_t fx_im_in = input_stage<In_W>(im_in[i], in_scale_factor);

        const out_t fx_norm = proc->process_sqr(fx_re_in, fx_im_in);

        results[i] = float(fx_norm) * ot_scale_factor;
    }

    constexpr auto scale =
        [](float value) {
            return (float) float(out_t(value / ot_scale_factor)) * ot_scale_factor;
        };
    constexpr float initial_out = float(decltype(proc)::element_type::initial_out * ot_scale_factor);

    double total_sig, total_qfx_noise, total_rel_error, total_sqr_error, total_abs_error, total_log_error;
    error_compute(double(norms_out[0]),
                  double(results[0] - initial_out),
                  total_sig, total_qfx_noise, total_rel_error,
                  total_sqr_error, total_abs_error, total_log_error);

    // Initialization Startup
    for (int64_t i = 1; i < int64_t(q - 1); i++) {
        INFO("input no " << i);
        REQUIRE_THAT(results[i] - results[i - 1],
                     Catch::Matchers::WithinAbs(scale(norms_out[i]) - scale(norms_out[i - 1]), ot_scale_factor));
        error_compute(double(norms_out[i]),
                      double(results[i] - initial_out),
                      total_sig, total_qfx_noise, total_rel_error,
                      total_sqr_error, total_abs_error, total_log_error);
    }

    REQUIRE_THAT(results[q - 1] - results[q - 2] + initial_out,
                 Catch::Matchers::WithinAbs(scale(norms_out[q - 1]) - scale(norms_out[q - 2]), ot_scale_factor));
    error_compute(double(norms_out[q - 1]),
                  double(results[q - 1]),
                  total_sig, total_qfx_noise, total_rel_error,
                  total_sqr_error, total_abs_error, total_log_error);

    for (int64_t i = q; i < int64_t(results.size()); i++) {
        INFO("input no " << i);
        REQUIRE_THAT(results[i] - results[i - 1],
                     Catch::Matchers::WithinAbs(scale(norms_out[i]) - scale(norms_out[i - 1]), 2 * ot_scale_factor));
        error_compute(double(norms_out[i]),
                      double(results[i]),
                      total_sig, total_qfx_noise, total_rel_error,
                      total_sqr_error, total_abs_error, total_log_error);
    }

    const double mre   = total_rel_error / norms_out.size();
    const double mae   = total_abs_error / norms_out.size();
    const double rmsle = std::sqrt(total_log_error / norms_out.size());
    const double rmse  = std::sqrt(total_sqr_error / norms_out.size());
    const double snr   = 20 * std::log10(total_sig / total_qfx_noise);
    Catch::cout() << "  -- Mean Relative Error = " << mre << "\n  -- Mean Absolute Error = " << mae << "\n  -- RMSE = " << rmse
                  << "\n  -- RMSLE = " << rmsle << "\n  -- SNR (dB) = " << snr << std::endl;
}
