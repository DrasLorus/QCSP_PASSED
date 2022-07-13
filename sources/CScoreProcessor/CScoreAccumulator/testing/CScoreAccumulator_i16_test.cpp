#include "CScoreProcessor/CScoreAccumulator/CScoreAccumulator.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <matioCpp/matioCpp.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CScoreAccumulator int16_t works for high snr inputs (q: 64, N: 60)", "[scoreaccu][high][fixed]") {

    constexpr unsigned q     = 64;
    constexpr unsigned N     = 60;
    constexpr unsigned In_W  = 24; // High SNR => High signal power !
    constexpr unsigned In_I  = 20; // High SNR => High signal power !
    constexpr unsigned Out_W = (In_W + pow2_log2<q>());
    constexpr unsigned Out_I = (In_I + pow2_log2<q>());

    constexpr float in_scale_factor = float(1U << (In_W - In_I));
    constexpr float ot_scale_factor = 1. / float(1 << (Out_W - Out_I));

    matioCpp::File data_file("../data/test_data.mat", matioCpp::FileMode::ReadOnly);

    const matioCpp::Vector<float> cabs_in = data_file.read("cabs_max_l2_infdB_w1_q64_N60_0_n10").asVector<float>();
    if (cabs_in.size() == 0 || not cabs_in.isValid()) {
        throw "cabs_max_l2_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const matioCpp::Vector<float> score_out = data_file.read("score_l2_infdB_w1_q64_N60_0_n10").asVector<float>();
    if (score_out.size() == 0 || not score_out.isValid()) {
        throw "score_l2_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    CScoreAccumulator<N, q, int16_t> * proc = new CScoreAccumulator<N, q, int16_t>();

    vector<float> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const uint32_t fx_in = uint32_t(cabs_in[i] * in_scale_factor);

        const uint32_t fx_score = proc->process(fx_in);

        results[i] = float(double(fx_score) * ot_scale_factor);
    }

    constexpr float margin = 480.f * 1e-3f; // Max of this test vector is 480

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(score_out[i], margin));
    }

    delete proc;
}

TEST_CASE("CScoreAccumulator int16_t works for low snr inputs (q: 64, N: 60)", "[scoreaccu][low][fixed]") {

    constexpr unsigned q     = 64;
    constexpr unsigned N     = 60;
    constexpr unsigned In_W  = 24;
    constexpr unsigned In_I  = 8;
    constexpr unsigned Out_W = (In_W + pow2_log2<q>() - 1);
    constexpr unsigned Out_I = (In_I + pow2_log2<q>() - 1);

    constexpr float  in_scale_factor = float(1U << (In_W - In_I));
    constexpr double ot_scale_factor = 1. / double(1 << (Out_W - Out_I));

    matioCpp::File data_file("../data/test_data.mat", matioCpp::FileMode::ReadOnly);

    const matioCpp::Vector<float> cabs_in = data_file.read("cabs_max_l2_m10dB_w1_q64_N60_0_n30").asVector<float>();
    if (cabs_in.size() == 0 || not cabs_in.isValid()) {
        throw "cabs_max_l2_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const matioCpp::Vector<float> score_out = data_file.read("score_l2_m10dB_w1_q64_N60_0_n30").asVector<float>();
    if (score_out.size() == 0 || not score_out.isValid()) {
        throw "score_l2_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    CScoreAccumulator<N, q, int16_t> * proc = new CScoreAccumulator<N, q, int16_t>();

    vector<float> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const uint32_t fx_in = uint32_t(cabs_in[i] * in_scale_factor);

        const uint32_t fx_score = proc->process(fx_in);

        results[i] = float(double(fx_score) * ot_scale_factor);
    }

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(score_out[i], 1e-4f));
    }

    delete proc;
}
