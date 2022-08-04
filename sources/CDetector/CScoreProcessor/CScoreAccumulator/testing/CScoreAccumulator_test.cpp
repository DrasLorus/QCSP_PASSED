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

TEST_CASE("CScoreAccumulator works for high snr inputs (q: 64, N: 60)", "[scoreaccu][high]") {

    constexpr unsigned q = 64;
    constexpr unsigned N = 60;

    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);

    const vector<float> cabs_in = load_data_vector<float, float>(data_file, "cabs_max_l2_infdB_w1_q64_N60_0_n10");

    const vector<float> score_out = load_data_vector<float, float>(data_file, "score_l2_infdB_w1_q64_N60_0_n10");

    Mat_Close(data_file);

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

TEMPLATE_TEST_CASE("CScoreAccumulator works for low snr inputs (q: 64, N: 60)", "[scoreaccu][low][q64][N60]", float, double) {

    constexpr unsigned q = 64;
    constexpr unsigned N = 60;

    mat_t * data_file = Mat_Open("../data/test_data_w1_nofreq.mat", MAT_ACC_RDONLY);

    const vector<TestType> cabs_in = load_data_vector<TestType, float>(data_file, "cabs_max_l2_m10dB_w1_q64_N60_0_n30");

    const vector<TestType> score_out = load_data_vector<TestType, float>(data_file, "score_l2_m10dB_w1_q64_N60_0_n30");

    Mat_Close(data_file);

    CScoreAccumulator<N, q, TestType> * proc = new CScoreAccumulator<N, q, TestType>();

    vector<TestType> results(score_out.size(), 0.f);

    for (int64_t i = 0; i < int64_t(results.size()); i++) {
        results[i] = proc->process(cabs_in[i]);
    }

    for (int64_t i = 0; i < q * N; i++) {
        INFO("input no " << i);
        REQUIRE_THAT(results[i], Catch::Matchers::WithinAbs(score_out[i], TestType(1e-4f)));
    }

    for (int64_t i = q * N; i < int64_t(results.size()); i++) {
        INFO("input no " << i);
        REQUIRE_THAT(results[i], Catch::Matchers::WithinRel(score_out[i], TestType(1e-4f)));
    }

    delete proc;
}
