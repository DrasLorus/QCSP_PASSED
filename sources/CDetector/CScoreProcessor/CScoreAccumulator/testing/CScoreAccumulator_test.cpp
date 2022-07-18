#include "CScoreProcessor/CScoreAccumulator/CScoreAccumulator.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <matioCpp/matioCpp.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEST_CASE("CScoreAccumulator works for high snr inputs (q: 64, N: 60)", "[scoreaccu][high]") {

    constexpr unsigned q = 64;
    constexpr unsigned N = 60;

    matioCpp::File data_file("../data/test_data.mat", matioCpp::FileMode::ReadOnly);

    const matioCpp::Vector<float> cabs_in = data_file.read("cabs_max_l2_infdB_w1_q64_N60_0_n10").asVector<float>();
    if (cabs_in.size() == 0 || not cabs_in.isValid()) {
        throw "cabs_max_l2_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    const matioCpp::Vector<float> score_out = data_file.read("score_l2_infdB_w1_q64_N60_0_n10").asVector<float>();
    if (score_out.size() == 0 || not score_out.isValid()) {
        throw "score_l2_infdB_w1_q64_N60_0_n10 can't be loaded.";
    }

    CScoreAccumulator<N, q> * proc = new CScoreAccumulator<N, q>();

    vector<float> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        results[i] = proc->process(cabs_in[i]);
    }

    constexpr float margin = 480.f * 1e-4;

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(score_out[i], margin));
    }

    delete proc;
}

TEST_CASE("CScoreAccumulator works for low snr inputs (q: 64, N: 60)", "[scoreaccu][low]") {

    constexpr unsigned q = 64;
    constexpr unsigned N = 60;

    matioCpp::File data_file("../data/test_data.mat", matioCpp::FileMode::ReadOnly);

    const matioCpp::Vector<float> cabs_in = data_file.read("cabs_max_l2_m10dB_w1_q64_N60_0_n30").asVector<float>();
    if (cabs_in.size() == 0 || not cabs_in.isValid()) {
        throw "cabs_max_l2_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    const matioCpp::Vector<float> score_out = data_file.read("score_l2_m10dB_w1_q64_N60_0_n30").asVector<float>();
    if (score_out.size() == 0 || not score_out.isValid()) {
        throw "score_l2_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    CScoreAccumulator<N, q> * proc = new CScoreAccumulator<N, q>();

    vector<float> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        results[i] = proc->process(cabs_in[i]);
    }

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(score_out[i], 1e-4f));
    }

    delete proc;
}
