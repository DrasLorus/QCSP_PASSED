#include "../CDetector.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <matio.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::vector;

TEMPLATE_TEST_CASE("CDetectorSerial (L2) works for high snr inputs (q: 64, N: 60, w: 1)", "[detector][high][l2]", float, double) {

    constexpr unsigned q       = 64;
    constexpr unsigned N       = 60;
    constexpr unsigned p_omega = 1;

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

    const vector<TestType> re_in((float *) cpx_ptr->Re, ((float *) cpx_ptr->Re) + tmp_mat->dims[0]);
    const vector<TestType> im_in((float *) cpx_ptr->Im, ((float *) cpx_ptr->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    // tmp_mat = Mat_VarRead(data_file, "score_l2_infdB_w1_q64_N60_0_n10");
    // if (not bool(tmp_mat)) {
    //     throw "score_l2_infdB_w1_q64_N60_0_n10 can't be loaded.";
    // }

    // const vector<float> score_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    // Mat_VarFree(tmp_mat);
    Mat_Close(data_file);
    constexpr TestType threshold        = .5f;
    constexpr uint32_t step_denominator = 1U;

    using detector_t = CDetectorSerial<N, q, p_omega, TestType, true>;

    detector_t * proc = new detector_t(
        pn.data(),
        threshold,
        step_denominator);

    REQUIRE(proc->threshold() == threshold);

    vector<std::pair<DetectionState<TestType, TestType, p_omega>, int64_t>> results;

    DetectionState<TestType, TestType, p_omega> running_state {false, false, {0.f}, 0.f, 0LU, 0LU, 0.f, 0U};

    SECTION("Standard Score") {
        for (int64_t i = 0; i < int64_t(re_in.size()); i++) {
            proc->process(re_in[i], im_in[i], &running_state);
            if (running_state.frame_detected and (running_state.chip_since_last_det == 0)) {
                results.push_back({running_state, i});
            }
            if (running_state.max_found and (running_state.chip_since_last_det == detector_t::window_size)) {
                results.push_back({running_state, i});
            }

            REQUIRE(running_state.max_found xor (running_state.chip_since_last_det < detector_t::window_size));
            REQUIRE(running_state.chip_since_last_det < (detector_t::window_size * 2));
            REQUIRE((running_state.chip_from_max < detector_t::window_size));
            REQUIRE(running_state.max_score > 0.f - 4.8e-4);
            REQUIRE((running_state.max_score == running_state.scores[0] or (running_state.frame_detected)));
        }

        REQUIRE(results.size() == 20LU);

        constexpr TestType low_score  = TestType(1.);
        constexpr TestType high_score = TestType(480.);
        constexpr TestType cfos       = TestType(0.);
        constexpr TestType tolerance  = TestType(1e-5);

        for (int64_t i = 0; i < int64_t(results.size()); i += 2) {
            const auto & result_i   = results[i];
            const auto & result_ip1 = results[i + 1];

            const int64_t chip_no = result_ip1.second - detector_t::window_size + result_ip1.first.chip_from_max;

            // * If you want to explore the results, uncomment the 6 following lines
            // printf("-- Detection no %ld score exceeded the threshold of %5.1f at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
            //        i >> 1, proc->threshold(), result_i.second, result_i.first.max_score, result_i.first.frequency_offset);
            // printf("-- Detection no %ld max has been found at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
            //        i >> 1, result_ip1.second - result_ip1.first.chip_from_max, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE(result_i.second == (i >> 1) * 19200 + 7680);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(low_score, tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(chip_no == ((i >> 1) * 19200 + 11520 - 1));
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(high_score, tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
    }

    SECTION("Square Score") {
        for (int64_t i = 0; i < int64_t(re_in.size()); i++) {
            proc->process_sqr(re_in[i], im_in[i], &running_state);
            if (running_state.frame_detected and (running_state.chip_since_last_det == 0)) {
                results.push_back({running_state, i});
            }
            if (running_state.max_found and (running_state.chip_since_last_det == detector_t::window_size)) {
                results.push_back({running_state, i});
            }

            REQUIRE(running_state.max_found xor (running_state.chip_since_last_det < detector_t::window_size));
            REQUIRE(running_state.chip_since_last_det < (detector_t::window_size * 2));
            REQUIRE((running_state.chip_from_max < detector_t::window_size));
            REQUIRE(running_state.max_score > 0.f - 4.8e-4);
            REQUIRE((running_state.max_score == running_state.scores[0] or (running_state.frame_detected)));
        }

        REQUIRE(results.size() == 20LU);

        constexpr TestType low_score  = TestType(1.);
        constexpr TestType high_score = TestType(3840.);
        constexpr TestType cfos       = TestType(0.);
        constexpr TestType tolerance  = TestType(1e-5);

        for (int64_t i = 0; i < int64_t(results.size()); i += 2) {
            const auto & result_i   = results[i];
            const auto & result_ip1 = results[i + 1];

            const int64_t chip_no = result_ip1.second - detector_t::window_size + result_ip1.first.chip_from_max;

            // * If you want to explore the results, uncomment the 6 following lines
            // printf("-- Detection no %ld score exceeded the threshold of %5.1f at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
            //        i >> 1, proc->threshold(), result_i.second, result_i.first.max_score, result_i.first.frequency_offset);
            // printf("-- Detection no %ld max has been found at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
            //        i >> 1, result_ip1.second - result_ip1.first.chip_from_max, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE(result_i.second == (i >> 1) * 19200 + 7680);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(low_score, tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(chip_no == ((i >> 1) * 19200 + 11520 - 1));
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(high_score, tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
    }
    delete proc;
}

TEMPLATE_TEST_CASE("CDetectorSerial (raw) works for high snr inputs (q: 64, N: 60, w: 1)", "[detector][high][raw]", float, double) {

    constexpr unsigned q       = 64;
    constexpr unsigned N       = 60;
    constexpr unsigned p_omega = 1;

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

    const vector<TestType> re_in((float *) cpx_ptr->Re, ((float *) cpx_ptr->Re) + tmp_mat->dims[0]);
    const vector<TestType> im_in((float *) cpx_ptr->Im, ((float *) cpx_ptr->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    // tmp_mat = Mat_VarRead(data_file, "score_l2_infdB_w1_q64_N60_0_n10");
    // if (not bool(tmp_mat)) {
    //     throw "score_l2_infdB_w1_q64_N60_0_n10 can't be loaded.";
    // }

    // const vector<float> score_out((float *) tmp_mat->data, ((float *) tmp_mat->data) + tmp_mat->dims[0]);

    // Mat_VarFree(tmp_mat);
    Mat_Close(data_file);
    constexpr TestType threshold        = .5f;
    constexpr uint32_t step_denominator = 1U;

    using detector_t = CDetectorSerial<N, q, p_omega, TestType, false>;

    detector_t * proc = new detector_t(
        pn.data(),
        threshold,
        step_denominator);

    REQUIRE(proc->threshold() == threshold);

    vector<std::pair<DetectionState<TestType, TestType, p_omega>, int64_t>> results;

    DetectionState<TestType, TestType, p_omega> running_state {false, false, {0.f}, 0.f, 0LU, 0LU, 0.f, 0U};

    SECTION("Standard Score") {
        for (int64_t i = 0; i < int64_t(re_in.size()); i++) {
            proc->process(re_in[i], im_in[i], &running_state);
            if (running_state.frame_detected and (running_state.chip_since_last_det == 0)) {
                results.push_back({running_state, i});
            }
            if (running_state.max_found and (running_state.chip_since_last_det == detector_t::window_size)) {
                results.push_back({running_state, i});
            }

            REQUIRE(running_state.max_found xor (running_state.chip_since_last_det < detector_t::window_size));
            REQUIRE(running_state.chip_since_last_det < (detector_t::window_size * 2));
            REQUIRE((running_state.chip_from_max < detector_t::window_size));
            REQUIRE(running_state.max_score > 0.f - 4.8e-4);
            REQUIRE((running_state.max_score == running_state.scores[0] or (running_state.frame_detected)));
        }

        REQUIRE(results.size() == 20LU);

        constexpr TestType low_score  = TestType(1.);
        constexpr TestType high_score = TestType(3840.);
        constexpr TestType cfos       = TestType(0.);
        constexpr TestType tolerance  = TestType(1e-5);

        for (int64_t i = 0; i < int64_t(results.size()); i += 2) {
            const auto & result_i   = results[i];
            const auto & result_ip1 = results[i + 1];

            const int64_t chip_no = result_ip1.second - detector_t::window_size + result_ip1.first.chip_from_max;

            // * If you want to explore the results, uncomment the 6 following lines
            // printf("-- Detection no %ld score exceeded the threshold of %5.1f at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
            //        i >> 1, proc->threshold(), result_i.second, result_i.first.max_score, result_i.first.frequency_offset);
            // printf("-- Detection no %ld max has been found at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
            //        i >> 1, result_ip1.second - result_ip1.first.chip_from_max, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE(result_i.second == (i >> 1) * 19200 + 7680);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(low_score, tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(chip_no == ((i >> 1) * 19200 + 11520 - 1));
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(high_score, tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
    }

    SECTION("Square Score") {
        for (int64_t i = 0; i < int64_t(re_in.size()); i++) {
            proc->process_sqr(re_in[i], im_in[i], &running_state);
            if (running_state.frame_detected and (running_state.chip_since_last_det == 0)) {
                results.push_back({running_state, i});
            }
            if (running_state.max_found and (running_state.chip_since_last_det == detector_t::window_size)) {
                results.push_back({running_state, i});
            }

            REQUIRE(running_state.max_found xor (running_state.chip_since_last_det < detector_t::window_size));
            REQUIRE(running_state.chip_since_last_det < (detector_t::window_size * 2));
            REQUIRE((running_state.chip_from_max < detector_t::window_size));
            REQUIRE(running_state.max_score > 0.f - 4.8e-4);
            REQUIRE((running_state.max_score == running_state.scores[0] or (running_state.frame_detected)));
        }

        REQUIRE(results.size() == 20LU);

        constexpr TestType low_score  = TestType(1.);
        constexpr TestType high_score = TestType(245760.);
        constexpr TestType cfos       = TestType(0.);
        constexpr TestType tolerance  = TestType(1e-5);

        for (int64_t i = 0; i < int64_t(results.size()); i += 2) {
            const auto & result_i   = results[i];
            const auto & result_ip1 = results[i + 1];

            const int64_t chip_no = result_ip1.second - detector_t::window_size + result_ip1.first.chip_from_max;

            // * If you want to explore the results, uncomment the 6 following lines
            // printf("-- Detection no %ld score exceeded the threshold of %5.1f at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
            //        i >> 1, proc->threshold(), result_i.second, result_i.first.max_score, result_i.first.frequency_offset);
            // printf("-- Detection no %ld max has been found at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
            //        i >> 1, result_ip1.second - result_ip1.first.chip_from_max, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE(result_i.second == (i >> 1) * 19200 + 7680);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(low_score, tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(chip_no == ((i >> 1) * 19200 + 11520 - 1));
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(high_score, tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
    }
    delete proc;
}

#if 0
TEST_CASE("CDetectorSerial works for low snr inputs (q: 64, N: 60)", "[scoreproc][low]") {

    constexpr unsigned q = 64;
    constexpr unsigned N = 60;

    matioCpp::File parm_file("../data/parameters_20210903.mat", MAT_ACC_RDONLY);

    matioCpp::Vector<double> PN64 = parm_file.read("PN64").asVector<double>();
    vector<float>            pn(PN64.size());

    std::copy(PN64.begin(), PN64.end(), pn.begin());

    if (PN64.size() != q) {
        throw "PN64 and q don't match.";
    }

    if (unsigned(parm_file.read("n_frame").asElement<double>()) != N) {
        throw "Fetched N and local N don't match.";
    }

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

    CDetectorSerial<N, q, float> * proc = new CDetectorSerial<N, q, float>(pn.data());

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
#endif
