#include "../CScoreAccumulator.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <matio.h>
#include <string>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::string;
using std::vector;

namespace {
template <typename T, typename data_type>
vector<T> load_data_vector(mat_t * matfile, const string & varname) {
    matvar_t * tmp_var = Mat_VarRead(matfile, varname.c_str());
    if (not bool(tmp_var)) {
        throw(varname + " can't be loaded.");
    }

    const vector<T> vect((data_type *) tmp_var->data, (data_type *) tmp_var->data + tmp_var->dims[0]);

    Mat_VarFree(tmp_var);

    return vect;
}
} // namespace

TEST_CASE("CScoreAccumulator int16_t works for high snr inputs (q: 64, N: 60)", "[scoreaccu][high][fixed]") {

    constexpr unsigned q     = 64;
    constexpr unsigned N     = 60;
    constexpr unsigned In_W  = 24; // High SNR => High signal power !
    constexpr unsigned In_I  = 20; // High SNR => High signal power !
    constexpr unsigned Out_W = (In_W + pow2_log2<q>());
    constexpr unsigned Out_I = (In_I + pow2_log2<q>());

    constexpr float in_scale_factor = float(1U << (In_W - In_I));
    constexpr float ot_scale_factor = 1. / float(1 << (Out_W - Out_I));

    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);

    const vector<float> cabs_in = load_data_vector<float, float>(data_file, "cabs_max_l2_infdB_w1_q64_N60_0_n10");

    const vector<float> score_out = load_data_vector<float, float>(data_file, "score_l2_infdB_w1_q64_N60_0_n10");

    Mat_Close(data_file);

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

    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);

    const vector<float> cabs_in = load_data_vector<float, float>(data_file, "cabs_max_l2_infdB_w1_q64_N60_0_n10");

    const vector<float> score_out = load_data_vector<float, float>(data_file, "score_l2_infdB_w1_q64_N60_0_n10");

    Mat_Close(data_file);

    CScoreAccumulator<N, q, int16_t> * proc = new CScoreAccumulator<N, q, int16_t>();

    vector<float> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        const uint32_t fx_in = uint32_t(cabs_in[i] * in_scale_factor);

        const uint32_t fx_score = proc->process(fx_in);

        results[i] = float(double(fx_score) * ot_scale_factor);
    }

    for (int64_t i = q * N; i < int64_t(results.size()); i++) {
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(score_out[i], 5e-4f));
    }

    delete proc;
}
