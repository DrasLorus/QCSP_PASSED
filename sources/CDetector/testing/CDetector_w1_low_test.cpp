#include "CDetector/CDetector.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <matio.h>
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

    const vector<T> vect((data_type *) tmp_var->data, ((data_type *) tmp_var->data) + tmp_var->dims[0]);

    Mat_VarFree(tmp_var);

    return vect;
}
} // namespace

// Doesn't work for double, as the score calculations differs.

TEST_CASE("CDetectorSerial (L2, float) works for low snr inputs (q: 64, N: 60, p_omega: 1)", "[detector][low][l2]") {

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

    matvar_t * tmp_mat = Mat_VarRead(data_file, "data_input_m10dB_w1_q64_N60_0_n30");
    if (not bool(tmp_mat)) {
        throw "data_input_m10dB_w1_q64_N60_0_n30 can't be loaded.";
    }

    mat_complex_split_t * cpx_ptr = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) cpx_ptr->Re, ((float *) cpx_ptr->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) cpx_ptr->Im, ((float *) cpx_ptr->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    constexpr unsigned run_size = N * q * 5;

    constexpr uint32_t step_denominator = 1U;

    using detector_t = CDetectorSerial<N, q, p_omega, float, true>;

    vector<std::pair<DetectionState<float, float, p_omega>, int64_t>> results;

    DetectionState<float, float, p_omega> running_state {false, false, {0.f}, 0.f, 0LU, 0LU, 0.f};

    SECTION("Standard Score") {
        constexpr float threshold = 140.5f;

        const vector<float> score_out = load_data_vector<float, float>(data_file, "score_l2_m10dB_w1_q64_N60_0_n30");

        vector<float>   det_scores(30, 0);
        vector<int64_t> det_idxs(30, 0);
        vector<float>   max_scores(30, 0);
        vector<int64_t> max_idxs(30, 0);
        for (unsigned u = 0; u < 30U; u++) {
            const auto begin = score_out.data() + u * run_size;
            const auto end   = score_out.data() + (u + 1) * run_size;

            const auto det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
            det_scores[u]      = *det_ptr;
            det_idxs[u]        = det_ptr - score_out.data();

            const auto max_ptr = std::max_element(begin, end);
            max_scores[u]      = *max_ptr;
            max_idxs[u]        = max_ptr - score_out.data();
        }
        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        REQUIRE(proc->threshold() == threshold);

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

        REQUIRE(results.size() == 60LU);

        constexpr float cfos      = float(0.);
        constexpr float tolerance = float(5e-5);

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
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE(result_i.second == det_idxs[i >> 1]);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(det_scores[i >> 1], tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(chip_no == max_idxs[(i >> 1)]);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(max_scores[(i >> 1)], tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
        delete proc;
    }

    SECTION("Square Score") {
        constexpr float threshold = float(350.5);

        const vector<float> score_out = load_data_vector<float, float>(data_file, "score_sqr_l2_m10dB_w1_q64_N60_0_n30");

        vector<float>   det_scores(30, 0);
        vector<int64_t> det_idxs(30, 0);
        vector<float>   max_scores(30, 0);
        vector<int64_t> max_idxs(30, 0);
        for (unsigned u = 0; u < 30U; u++) {
            const auto begin = score_out.data() + u * run_size;
            const auto end   = score_out.data() + (u + 1) * run_size;

            const auto det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
            det_scores[u]      = *det_ptr;
            det_idxs[u]        = det_ptr - score_out.data();

            const auto max_ptr = std::max_element(begin, end);
            max_scores[u]      = *max_ptr;
            max_idxs[u]        = max_ptr - score_out.data();
        }
        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        REQUIRE(proc->threshold() == threshold);

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

        REQUIRE(results.size() == 60LU);

        constexpr float cfos      = float(0.);
        constexpr float tolerance = float(5e-5);

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
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE((result_i.second == det_idxs[i >> 1]));
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(det_scores[i >> 1], tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(chip_no == max_idxs[(i >> 1)]);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(max_scores[(i >> 1)], tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
        delete proc;
    }
    Mat_Close(data_file);
}

TEST_CASE("CDetectorSerial (L2, double) works for low snr inputs (q: 64, N: 60, p_omega: 1)", "[detector][low][l2]") {

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

    const vector<double> pn((double *) tmp_pn->data, (double *) tmp_pn->data + tmp_pn->dims[1]);
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

    const vector<double> re_in((float *) cpx_ptr->Re, ((float *) cpx_ptr->Re) + tmp_mat->dims[0]);
    const vector<double> im_in((float *) cpx_ptr->Im, ((float *) cpx_ptr->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);
    constexpr unsigned run_size = N * q * 5;

    constexpr uint32_t step_denominator = 1U;

    using detector_t = CDetectorSerial<N, q, p_omega, double, true>;

    vector<std::pair<DetectionState<double, double, p_omega>, int64_t>> results;

    DetectionState<double, double, p_omega> running_state {false, false, {0.}, 0., 0LU, 0LU, 0.};

    SECTION("Standard Score") {
        constexpr double threshold = 140.5f;

        const vector<double> score_out = load_data_vector<double, float>(data_file, "score_l2_m10dB_w1_q64_N60_0_n30");

        vector<int64_t> det_idxs(30, 0);
        vector<int64_t> max_idxs(30, 0);
        for (unsigned u = 0; u < 30U; u++) {
            const double * begin = score_out.data() + u * run_size;
            const double * end   = score_out.data() + (u + 1) * run_size;

            const double * det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
            det_idxs[u]            = det_ptr - score_out.data();

            const double * max_ptr = std::max_element(begin, end);
            max_idxs[u]            = max_ptr - score_out.data();
        }
        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        REQUIRE(proc->threshold() == threshold);

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

        REQUIRE(results.size() == 60LU);

        constexpr double cfos      = double(0.);
        constexpr double tolerance = double(1e-5);

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
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE((std::abs(result_i.second - det_idxs[i >> 1]) < 10));
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(score_out[result_i.second], tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(std::abs(chip_no - max_idxs[(i >> 1)]) < 10);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(score_out[chip_no], tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
        delete proc;
    }

    SECTION("Square Score") {
        constexpr double threshold = double(350.5);

        const vector<double> score_out = load_data_vector<double, float>(data_file, "score_sqr_l2_m10dB_w1_q64_N60_0_n30");

        vector<int64_t> det_idxs(30, 0);
        vector<int64_t> max_idxs(30, 0);
        for (unsigned u = 0; u < 30U; u++) {
            const auto begin = score_out.data() + u * run_size;
            const auto end   = score_out.data() + (u + 1) * run_size;

            const auto det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
            det_idxs[u]        = det_ptr - score_out.data();

            const auto max_ptr = std::max_element(begin, end);
            max_idxs[u]        = max_ptr - score_out.data();
        }
        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        REQUIRE(proc->threshold() == threshold);

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

        REQUIRE(results.size() == 60LU);

        constexpr double cfos      = double(0.);
        constexpr double tolerance = double(5e-5);

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
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE((std::abs(result_i.second - det_idxs[i >> 1]) < 10));
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(score_out[result_i.second], tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(std::abs(chip_no - max_idxs[(i >> 1)]) < 10);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(score_out[chip_no], tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
        delete proc;
    }
    Mat_Close(data_file);
}

TEST_CASE("CDetectorSerial (Raw, float) works for low snr inputs (q: 64, N: 60, p_omega: 1)", "[detector][low][raw]") {

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

    constexpr unsigned run_size = N * q * 5;

    constexpr uint32_t step_denominator = 1U;

    using detector_t = CDetectorSerial<N, q, p_omega, float, false>;

    vector<std::pair<DetectionState<float, float, p_omega>, int64_t>> results;

    DetectionState<float, float, p_omega> running_state {false, false, {0.f}, 0.f, 0LU, 0LU, 0.f};

    SECTION("Standard Score") {
        constexpr float threshold = 3600.5f;

        const vector<float> score_out = load_data_vector<float, float>(data_file, "score_raw_m10dB_w1_q64_N60_0_n30");

        vector<float>   det_scores(30, 0);
        vector<int64_t> det_idxs(30, 0);
        vector<float>   max_scores(30, 0);
        vector<int64_t> max_idxs(30, 0);
        for (unsigned u = 0; u < 30U; u++) {
            const auto begin = score_out.data() + u * run_size;
            const auto end   = score_out.data() + (u + 1) * run_size;

            const auto det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
            det_scores[u]      = *det_ptr;
            det_idxs[u]        = det_ptr - score_out.data();

            const auto max_ptr = std::max_element(begin, end);
            max_scores[u]      = *max_ptr;
            max_idxs[u]        = max_ptr - score_out.data();
        }
        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        REQUIRE(proc->threshold() == threshold);

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

        REQUIRE(results.size() == 60LU);

        constexpr float cfos      = float(0.);
        constexpr float tolerance = float(5e-5);

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
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE(result_i.second == det_idxs[i >> 1]);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(det_scores[i >> 1], tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(chip_no == max_idxs[(i >> 1)]);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(max_scores[(i >> 1)], tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
        delete proc;
    }

    SECTION("Square Score") {
        constexpr float threshold = 220500.f;

        const vector<float> score_out = load_data_vector<float, float>(data_file, "score_sqr_raw_m10dB_w1_q64_N60_0_n30");

        vector<float>   det_scores(30, 0);
        vector<int64_t> det_idxs(30, 0);
        vector<float>   max_scores(30, 0);
        vector<int64_t> max_idxs(30, 0);
        for (unsigned u = 0; u < 30U; u++) {
            const auto begin = score_out.data() + u * run_size;
            const auto end   = score_out.data() + (u + 1) * run_size;

            const auto det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
            det_scores[u]      = *det_ptr;
            det_idxs[u]        = det_ptr - score_out.data();

            const auto max_ptr = std::max_element(begin, end);
            max_scores[u]      = *max_ptr;
            max_idxs[u]        = max_ptr - score_out.data();
        }
        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        REQUIRE(proc->threshold() == threshold);

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

        REQUIRE(results.size() == 60LU);

        constexpr float cfos      = float(0.);
        constexpr float tolerance = float(5e-5);

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
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE((result_i.second == det_idxs[i >> 1]));
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(det_scores[i >> 1], tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(chip_no == max_idxs[(i >> 1)]);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(max_scores[(i >> 1)], tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
        delete proc;
    }
    Mat_Close(data_file);
}

TEST_CASE("CDetectorSerial (Raw, double) works for low snr inputs (q: 64, N: 60, p_omega: 1)", "[detector][low][raw]") {

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

    const vector<double> pn((double *) tmp_pn->data, (double *) tmp_pn->data + tmp_pn->dims[1]);
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

    const vector<double> re_in((float *) cpx_ptr->Re, ((float *) cpx_ptr->Re) + tmp_mat->dims[0]);
    const vector<double> im_in((float *) cpx_ptr->Im, ((float *) cpx_ptr->Im) + tmp_mat->dims[0]);

    Mat_VarFree(tmp_mat);

    constexpr unsigned run_size = N * q * 5;

    constexpr uint32_t step_denominator = 1U;

    using detector_t = CDetectorSerial<N, q, p_omega, double, false>;

    vector<std::pair<DetectionState<double, double, p_omega>, int64_t>> results;

    DetectionState<double, double, p_omega> running_state {false, false, {0.}, 0., 0LU, 0LU, 0.};

    SECTION("Standard Score") {
        constexpr double threshold = 3600.5;

        const vector<double> score_out = load_data_vector<double, float>(data_file, "score_raw_m10dB_w1_q64_N60_0_n30");

        vector<int64_t> det_idxs(30, 0);
        vector<int64_t> max_idxs(30, 0);
        for (unsigned u = 0; u < 30U; u++) {
            const double * begin = score_out.data() + u * run_size;
            const double * end   = score_out.data() + (u + 1) * run_size;

            const double * det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
            det_idxs[u]            = det_ptr - score_out.data();

            const double * max_ptr = std::max_element(begin, end);
            max_idxs[u]            = max_ptr - score_out.data();
        }
        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        REQUIRE(proc->threshold() == threshold);

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

        REQUIRE(results.size() == 60LU);

        constexpr double cfos      = double(0.);
        constexpr double tolerance = double(1e-5);

        for (int64_t i = 0; i < int64_t(results.size()); i += 2) {
            const auto & result_i   = results[i];
            const auto & result_ip1 = results[i + 1];

            const int64_t chip_no = result_ip1.second - detector_t::window_size + result_ip1.first.chip_from_max;

            // * If you want to explore the results, uncomment the 6 following lines
            printf("-- Detection no %ld score exceeded the threshold of %5.1f at chip no %ld,\n"
                   "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
                   i >> 1, proc->threshold(), result_i.second, result_i.first.max_score, result_i.first.frequency_offset);
            printf("-- Detection no %ld max has been found at chip no %ld,\n"
                   "-- for a score of %16.8f and a frequency offset of %16.8e Hz.\n",
                   i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE((std::abs(result_i.second - det_idxs[i >> 1]) < 10));
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(score_out[result_i.second], tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(std::abs(chip_no - max_idxs[(i >> 1)]) < 10);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(score_out[chip_no], tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
        delete proc;
    }

    SECTION("Square Score") {
        constexpr double threshold = 220500.;

        const vector<double> score_out = load_data_vector<double, float>(data_file, "score_sqr_raw_m10dB_w1_q64_N60_0_n30");

        vector<int64_t> det_idxs(30, 0);
        vector<int64_t> max_idxs(30, 0);
        for (unsigned u = 0; u < 30U; u++) {
            const auto begin = score_out.data() + u * run_size;
            const auto end   = score_out.data() + (u + 1) * run_size;

            const auto det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
            det_idxs[u]        = det_ptr - score_out.data();

            const auto max_ptr = std::max_element(begin, end);
            max_idxs[u]        = max_ptr - score_out.data();
        }
        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        REQUIRE(proc->threshold() == threshold);

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

        REQUIRE(results.size() == 60LU);

        constexpr double cfos      = double(0.);
        constexpr double tolerance = double(5e-5);

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
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset);

            REQUIRE((std::abs(result_i.second - det_idxs[i >> 1]) < 10));
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(score_out[result_i.second], tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
            REQUIRE(std::abs(chip_no - max_idxs[(i >> 1)]) < 10);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(score_out[chip_no], tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(cfos, tolerance));
        }
        delete proc;
    }
    Mat_Close(data_file);
}