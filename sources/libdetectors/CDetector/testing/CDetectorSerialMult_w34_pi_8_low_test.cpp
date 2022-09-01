#include "../CDetector.hpp"

#include <algorithm>
#include <catch2/catch.hpp>
#include <iomanip>
#include <iostream>
#include <matio.h>
#include <vector>

using namespace QCSP::StandaloneDetector;
using std::cerr; // NOLINT
using std::endl; // NOLINT
using std::string;
using std::vector;

namespace {
template <typename T, typename data_type, unsigned dim = 0>
vector<T> load_data_vector(mat_t * matfile, const string & varname) {
    matvar_t * tmp_var = Mat_VarRead(matfile, varname.c_str());
    if (not bool(tmp_var)) {
        throw(varname + " can't be loaded.");
    }

    const vector<T> vect((data_type *) tmp_var->data, ((data_type *) tmp_var->data) + tmp_var->dims[dim]);

    Mat_VarFree(tmp_var);

    return vect;
}
} // namespace

namespace {
template <typename T, typename data_type>
vector<vector<T>> load_data_matrix(mat_t * matfile, const string & varname) {
    matvar_t * tmp_var = Mat_VarRead(matfile, varname.c_str());
    if (not bool(tmp_var)) {
        throw(varname + " can't be loaded.");
    }
    vector<vector<T>> mat(tmp_var->dims[1]);
    for (size_t i = 0; i < tmp_var->dims[1]; i++) {
        mat[i] = vector<T>((data_type *) tmp_var->data + i * tmp_var->dims[0], ((data_type *) tmp_var->data) + (i + 1) * tmp_var->dims[0]);
    }

    Mat_VarFree(tmp_var);

    return mat;
}
} // namespace

TEST_CASE("CDetectorSerialMult (L2, float) works for low snr inputs (q: 64, N: 60, w: 34, d: 8)", "[detector][low][l2][w34][mult]") {

    constexpr unsigned q                = 64;
    constexpr unsigned N                = 60;
    constexpr unsigned p_omega          = 34;
    constexpr unsigned step_denominator = 8;

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

    mat_t * data_file = Mat_Open("../data/test_data_q64_N60_w34_step8_span2.0.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data_q64_N60_w34_step8_span2.0.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "data_input_m10dB_w34_q64_N60_2pi_1_n30");
    if (not bool(tmp_mat)) {
        throw "data_input_m10dB_w34_q64_N60_2pi_1_n30 can't be loaded.";
    }

    mat_complex_split_t * cpx_ptr = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) cpx_ptr->Re, ((float *) cpx_ptr->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) cpx_ptr->Im, ((float *) cpx_ptr->Im) + tmp_mat->dims[0]);

    // for (int i = 0; i < (int) re_in.size(); i++) {
    //     cerr << std::setw(16) << std::setprecision(8) << re_in[i] << " " << im_in[i] << endl;
    // }

    Mat_VarFree(tmp_mat);

    const vector<float> rotations = load_data_vector<float, float>(data_file, "rotations_m10dB_w34_q64_N60_2pi_1_n30");

    constexpr unsigned run_size = N * q * 5;

    using detector_t = CDetectorSerialMult<N, q, p_omega, float, true>;

    vector<std::pair<DetectionState<float, float, p_omega>, int64_t>> results;

    DetectionState<float, float, p_omega> running_state {false, false, {0.f}, 0.f, 0LU, 0LU, 0.f, 0U};

    SECTION("Standard Score") {
        constexpr float threshold = 140.5f;

        const vector<vector<float>> scores_out
            = load_data_matrix<float, float>(data_file, "score_sqrt_l2_m10dB_w34_q64_N60_2pi_1_n30");

        vector<std::pair<int64_t, int64_t>> det_idxs(30, {0, 0});
        vector<std::pair<int64_t, int64_t>> max_idxs(30, {0, 0});
        for (unsigned u = 0; u < 30U; u++) {
            vector<int64_t> local_det_idx(p_omega);
            vector<float>   local_det_score(p_omega);
            vector<int64_t> local_max_idx(p_omega);
            vector<float>   local_max_score(p_omega);
            for (unsigned w = 0; w < p_omega; w++) {
                const float * begin = scores_out[w].data() + u * run_size;
                const float * end   = scores_out[w].data() + (u + 1) * run_size;

                const float * det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
                local_det_score[w]    = *det_ptr;
                local_det_idx[w]      = det_ptr - scores_out[w].data();

                const float * max_ptr = std::max_element(begin, end);
                local_max_score[w]    = *max_ptr;
                local_max_idx[w]      = max_ptr - scores_out[w].data();
            }

            // * Search the minimum detection index, and if there is several
            // * frequencies that exceed the threshold at the same time, the one
            // * corresponding to the highest score is kept.
            int64_t current_freq = 0;
            for (unsigned w = 0; w < p_omega; w++) {
                const int current_idx = local_det_idx[current_freq];
                const int local_idx   = local_det_idx[w];

                if (local_idx < current_idx) {
                    current_freq = w;
                } else if (local_idx == current_idx) {
                    current_freq = local_det_score[w] > local_det_score[current_freq] ? w : current_freq;

                    // * NOTE: Uncomment following lines to monitor this case (two frequency hypothesis exceed the threshold at the same time)
                    // printf("%3u | %3lu -- local %10d: %16.8f current %2d: %16.8f\n\n",
                    //        u, current_freq,
                    //        local_idx, local_det_score[w],
                    //        current_idx, local_det_score[current_freq]);
                }
            }
            const unsigned det_freq = current_freq;

            det_idxs[u] = {det_freq, local_det_idx[det_freq]};

            const float *  max_freq_score = std::max_element(local_max_score.data(), local_max_score.data() + p_omega);
            const unsigned max_freq       = max_freq_score - local_max_score.data();

            max_idxs[u] = {max_freq, local_max_idx[max_freq]};
        }

        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        for (unsigned i = 0; i < p_omega; i++) {
            const float proc_freq = proc->frequency_error(i) * 2 * q;
            const float ref_freq  = float(i) / step_denominator - float(p_omega - 1) / float(2 * step_denominator);
            REQUIRE_THAT(proc_freq, Catch::Matchers::WithinRel(ref_freq, 1e-6f));
        }

        REQUIRE(proc->threshold() == threshold);

        // * NOTE: Uncomment the following line and the two next groups to log the score
        // vector<vector<float>> result_scores(scores_out[0].size(), vector<float>(p_omega, 0));

        for (int64_t i = 0; i < int64_t(re_in.size()); i++) {
            proc->process(re_in[i], im_in[i], &running_state);

            // * NOTE: Uncomment the following line, the next group and the previous group to log the score
            // memcpy(result_scores[i].data(), running_state.scores, p_omega * sizeof(float));

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
            REQUIRE((running_state.max_score == running_state.scores[running_state.frequency_index] or (running_state.frame_detected)));
        }

        // * NOTE: Uncomment the following lines and the two previous groups to log the score
        // FILE * f = fopen("scores.dat", "w");
        // for (int j = 0; j < (int) result_scores.size(); j++) {
        //     fprintf(f, "%16.8f %16.8f %16.8f %16.8f\n", result_scores[j][0], result_scores[j][1], result_scores[j][2], result_scores[j][3]);
        // }
        // fclose(f);

        REQUIRE(results.size() == 60LU);

        constexpr float tolerance = float(5e-5);

        for (int64_t i = 0; i < int64_t(results.size()); i += 2) {
            const auto & result_i   = results[i];
            const auto & result_ip1 = results[i + 1];

            const int64_t chip_no = result_ip1.second - detector_t::window_size + result_ip1.first.chip_from_max;

            const float det_cfos = proc->frequency_error(det_idxs[i >> 1].first);
            const float max_cfos = proc->frequency_error(max_idxs[i >> 1].first);

            // * NOTE: If you want to explore the results, uncomment the 6 following lines
            // printf("-- Detection no %ld score exceeded the threshold of %5.1f at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e (%3d / %u) Hz.\n",
            //        i >> 1, proc->threshold(), result_i.second, result_i.first.max_score, result_i.first.frequency_offset,
            //        result_i.first.frequency_index, p_omega);
            // printf("-- Detection no %ld max has been found at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e (%3d / %u) Hz.\n",
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset,
            //        result_i.first.frequency_index, p_omega);

            const float det_score = scores_out[det_idxs[i >> 1].first][det_idxs[i >> 1].second];
            const float max_score = scores_out[max_idxs[i >> 1].first][max_idxs[i >> 1].second];

            REQUIRE(result_i.first.frequency_index == det_idxs[i >> 1].first);
            REQUIRE(result_i.second == det_idxs[i >> 1].second);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(det_score, tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(det_cfos, tolerance));
            REQUIRE(result_ip1.first.frequency_index == max_idxs[i >> 1].first);
            REQUIRE(chip_no == max_idxs[(i >> 1)].second);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(max_score, tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(max_cfos, tolerance));
        }
        delete proc;
    }

    SECTION("Square Score") {

        constexpr float threshold = float(350.5);

        const vector<vector<float>> scores_out
            = load_data_matrix<float, float>(data_file, "score_sqr_l2_m10dB_w34_q64_N60_2pi_1_n30");

        vector<std::pair<int64_t, int64_t>> det_idxs(30, {0, 0});
        vector<std::pair<int64_t, int64_t>> max_idxs(30, {0, 0});
        for (unsigned u = 0; u < 30U; u++) {
            vector<int64_t> local_det_idx(p_omega);
            vector<float>   local_det_score(p_omega);
            vector<int64_t> local_max_idx(p_omega);
            vector<float>   local_max_score(p_omega);
            for (unsigned w = 0; w < p_omega; w++) {
                const float * begin = scores_out[w].data() + u * run_size;
                const float * end   = scores_out[w].data() + (u + 1) * run_size;

                const float * det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
                local_det_score[w]    = *det_ptr;
                local_det_idx[w]      = det_ptr - scores_out[w].data();

                const float * max_ptr = std::max_element(begin, end);
                local_max_score[w]    = *max_ptr;
                local_max_idx[w]      = max_ptr - scores_out[w].data();
            }

            // * Search the minimum detection index, and if there is several
            // * frequencies that exceed the threshold at the same time, the one
            // * corresponding to the highest score is kept.
            int64_t current_freq = 0;
            for (unsigned w = 0; w < p_omega; w++) {
                const int current_idx = local_det_idx[current_freq];
                const int local_idx   = local_det_idx[w];

                if (local_idx < current_idx) {
                    current_freq = w;
                } else if (local_idx == current_idx) {
                    current_freq = local_det_score[w] > local_det_score[current_freq] ? w : current_freq;

                    // * NOTE: Uncomment following lines to monitor this case (two frequency hypothesis exceed the threshold at the same time)
                    // printf("%3u | %3lu -- local %10d: %16.8f current %2d: %16.8f\n\n",
                    //        u, current_freq,
                    //        local_idx, local_det_score[w],
                    //        current_idx, local_det_score[current_freq]);
                }
            }
            const unsigned det_freq = current_freq;

            det_idxs[u] = {det_freq, local_det_idx[det_freq]};

            const float *  max_freq_score = std::max_element(local_max_score.data(), local_max_score.data() + p_omega);
            const unsigned max_freq       = max_freq_score - local_max_score.data();

            max_idxs[u] = {max_freq, local_max_idx[max_freq]};
        }

        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        for (unsigned i = 0; i < p_omega; i++) {
            const float proc_freq = proc->frequency_error(i) * 2 * q;
            const float ref_freq  = float(i) / step_denominator - float(p_omega - 1) / float(2 * step_denominator);
            REQUIRE_THAT(proc_freq, Catch::Matchers::WithinRel(ref_freq, 1e-6f));
        }

        REQUIRE(proc->threshold() == threshold);

        // * NOTE: Uncomment the following line and the two next groups to log the score
        // vector<vector<float>> result_scores(scores_out[0].size(), vector<float>(p_omega, 0));

        for (int64_t i = 0; i < int64_t(re_in.size()); i++) {
            proc->process_sqr(re_in[i], im_in[i], &running_state);

            // * NOTE: Uncomment the following line, the next group and the previous group to log the score
            // memcpy(result_scores[i].data(), running_state.scores, p_omega * sizeof(float));

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
            REQUIRE((running_state.max_score == running_state.scores[running_state.frequency_index] or (running_state.frame_detected)));
        }

        // * NOTE: Uncomment the following lines and the two previous groups to log the score
        // FILE * f = fopen("scores.dat", "w");
        // for (int j = 0; j < (int) result_scores.size(); j++) {
        //     fprintf(f, "%16.8f %16.8f %16.8f %16.8f\n", result_scores[j][0], result_scores[j][1], result_scores[j][2], result_scores[j][3]);
        // }
        // fclose(f);

        REQUIRE(results.size() == 60LU);

        constexpr float tolerance = float(5e-5);

        for (int64_t i = 0; i < int64_t(results.size()); i += 2) {
            const auto & result_i   = results[i];
            const auto & result_ip1 = results[i + 1];

            const int64_t chip_no = result_ip1.second - detector_t::window_size + result_ip1.first.chip_from_max;

            const float det_cfos = proc->frequency_error(det_idxs[i >> 1].first);
            const float max_cfos = proc->frequency_error(max_idxs[i >> 1].first);

            // * NOTE: If you want to explore the results, uncomment the 6 following lines
            // printf("-- Detection no %ld score exceeded the threshold of %5.1f at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e (%3d / %u) Hz.\n",
            //        i >> 1, proc->threshold(), result_i.second, result_i.first.max_score, result_i.first.frequency_offset,
            //        result_i.first.frequency_index, p_omega);
            // printf("-- Detection no %ld max has been found at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e (%3d / %u) Hz.\n",
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset,
            //        result_i.first.frequency_index, p_omega);

            const float det_score = scores_out[det_idxs[i >> 1].first][det_idxs[i >> 1].second];
            const float max_score = scores_out[max_idxs[i >> 1].first][max_idxs[i >> 1].second];

            REQUIRE(result_i.first.frequency_index == det_idxs[i >> 1].first);
            REQUIRE(result_i.second == det_idxs[i >> 1].second);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(det_score, tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(det_cfos, tolerance));
            REQUIRE(result_ip1.first.frequency_index == max_idxs[i >> 1].first);
            REQUIRE(chip_no == max_idxs[(i >> 1)].second);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(max_score, tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(max_cfos, tolerance));
        }
        delete proc;
    }
    Mat_Close(data_file);
}

TEST_CASE("CDetectorSerialMult (Raw, float) works for low snr inputs (q: 64, N: 60, w: 34, d: 8)", "[detector][low][raw][w34][mult]") {

    constexpr unsigned q                = 64;
    constexpr unsigned N                = 60;
    constexpr unsigned p_omega          = 34;
    constexpr unsigned step_denominator = 8;

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

    mat_t * data_file = Mat_Open("../data/test_data_q64_N60_w34_step8_span2.0.mat", MAT_ACC_RDONLY);
    if (not bool(data_file)) {
        throw "../data/test_data_q64_N60_w34_step8_span2.0.mat can't be opened.";
    }

    matvar_t * tmp_mat = Mat_VarRead(data_file, "data_input_m10dB_w34_q64_N60_2pi_1_n30");
    if (not bool(tmp_mat)) {
        throw "data_input_m10dB_w34_q64_N60_2pi_1_n30 can't be loaded.";
    }

    mat_complex_split_t * cpx_ptr = (mat_complex_split_t *) tmp_mat->data;

    const vector<float> re_in((float *) cpx_ptr->Re, ((float *) cpx_ptr->Re) + tmp_mat->dims[0]);
    const vector<float> im_in((float *) cpx_ptr->Im, ((float *) cpx_ptr->Im) + tmp_mat->dims[0]);

    // for (int i = 0; i < (int) re_in.size(); i++) {
    //     cerr << std::setw(16) << std::setprecision(8) << re_in[i] << " " << im_in[i] << endl;
    // }

    Mat_VarFree(tmp_mat);

    const vector<float> rotations = load_data_vector<float, float>(data_file, "rotations_m10dB_w34_q64_N60_2pi_1_n30");

    constexpr unsigned run_size = N * q * 5;

    using detector_t = CDetectorSerialMult<N, q, p_omega, float, false>;

    vector<std::pair<DetectionState<float, float, p_omega>, int64_t>> results;

    DetectionState<float, float, p_omega> running_state {false, false, {0.f}, 0.f, 0LU, 0LU, 0.f, 0U};

    SECTION("Standard Score") {
        constexpr float threshold = 3600.5f;

        const vector<vector<float>> scores_out
            = load_data_matrix<float, float>(data_file, "score_sqrt_raw_m10dB_w34_q64_N60_2pi_1_n30");

        vector<std::pair<int64_t, int64_t>> det_idxs(30, {0, 0});
        vector<std::pair<int64_t, int64_t>> max_idxs(30, {0, 0});
        for (unsigned u = 0; u < 30U; u++) {
            vector<int64_t> local_det_idx(p_omega);
            vector<float>   local_det_score(p_omega);
            vector<int64_t> local_max_idx(p_omega);
            vector<float>   local_max_score(p_omega);
            for (unsigned w = 0; w < p_omega; w++) {
                const float * begin = scores_out[w].data() + u * run_size;
                const float * end   = scores_out[w].data() + (u + 1) * run_size;

                const float * det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
                local_det_score[w]    = *det_ptr;
                local_det_idx[w]      = det_ptr - scores_out[w].data();

                const float * max_ptr = std::max_element(begin, end);
                local_max_score[w]    = *max_ptr;
                local_max_idx[w]      = max_ptr - scores_out[w].data();
            }

            // * Search the minimum detection index, and if there is several
            // * frequencies that exceed the threshold at the same time, the one
            // * corresponding to the highest score is kept.
            int64_t current_freq = 0;
            for (unsigned w = 0; w < p_omega; w++) {
                const int current_idx = local_det_idx[current_freq];
                const int local_idx   = local_det_idx[w];

                if (local_idx < current_idx) {
                    current_freq = w;
                } else if (local_idx == current_idx) {
                    current_freq = local_det_score[w] > local_det_score[current_freq] ? w : current_freq;

                    // * NOTE: Uncomment following lines to monitor this case (two frequency hypothesis exceed the threshold at the same time)
                    // printf("%3u | %3lu -- local %10d: %16.8f current %2d: %16.8f\n\n",
                    //        u, current_freq,
                    //        local_idx, local_det_score[w],
                    //        current_idx, local_det_score[current_freq]);
                }
            }
            const unsigned det_freq = current_freq;

            det_idxs[u] = {det_freq, local_det_idx[det_freq]};

            const float *  max_freq_score = std::max_element(local_max_score.data(), local_max_score.data() + p_omega);
            const unsigned max_freq       = max_freq_score - local_max_score.data();

            max_idxs[u] = {max_freq, local_max_idx[max_freq]};
        }

        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        for (unsigned i = 0; i < p_omega; i++) {
            const float proc_freq = proc->frequency_error(i) * 2 * q;
            const float ref_freq  = float(i) / step_denominator - float(p_omega - 1) / float(2 * step_denominator);
            REQUIRE_THAT(proc_freq, Catch::Matchers::WithinRel(ref_freq, 1e-6f));
        }

        REQUIRE(proc->threshold() == threshold);

        // * NOTE: Uncomment the following line and the two next groups to log the score
        // vector<vector<float>> result_scores(scores_out[0].size(), vector<float>(p_omega, 0));

        for (int64_t i = 0; i < int64_t(re_in.size()); i++) {
            proc->process(re_in[i], im_in[i], &running_state);

            // * NOTE: Uncomment the following line, the next group and the previous group to log the score
            // memcpy(result_scores[i].data(), running_state.scores, p_omega * sizeof(float));

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
            REQUIRE((running_state.max_score == running_state.scores[running_state.frequency_index] or (running_state.frame_detected)));
        }

        // * NOTE: Uncomment the following lines and the two previous groups to log the score
        // FILE * f = fopen("scores.dat", "w");
        // for (int j = 0; j < (int) result_scores.size(); j++) {
        //     fprintf(f, "%16.8f %16.8f %16.8f %16.8f\n", result_scores[j][0], result_scores[j][1], result_scores[j][2], result_scores[j][3]);
        // }
        // fclose(f);

        REQUIRE(results.size() == 60LU);

        constexpr float tolerance = float(5e-5);

        for (int64_t i = 0; i < int64_t(results.size()); i += 2) {
            const auto & result_i   = results[i];
            const auto & result_ip1 = results[i + 1];

            const int64_t chip_no = result_ip1.second - detector_t::window_size + result_ip1.first.chip_from_max;

            const float det_cfos = proc->frequency_error(det_idxs[i >> 1].first);
            const float max_cfos = proc->frequency_error(max_idxs[i >> 1].first);

            // * NOTE: If you want to explore the results, uncomment the 6 following lines
            // printf("-- Detection no %ld score exceeded the threshold of %5.1f at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e (%3d / %u) Hz.\n",
            //        i >> 1, proc->threshold(), result_i.second, result_i.first.max_score, result_i.first.frequency_offset,
            //        result_i.first.frequency_index, p_omega);
            // printf("-- Detection no %ld max has been found at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e (%3d / %u) Hz.\n",
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset,
            //        result_i.first.frequency_index, p_omega);

            const float det_score = scores_out[det_idxs[i >> 1].first][det_idxs[i >> 1].second];
            const float max_score = scores_out[max_idxs[i >> 1].first][max_idxs[i >> 1].second];

            REQUIRE(result_i.first.frequency_index == det_idxs[i >> 1].first);
            REQUIRE(result_i.second == det_idxs[i >> 1].second);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(det_score, tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(det_cfos, tolerance));
            REQUIRE(result_ip1.first.frequency_index == max_idxs[i >> 1].first);
            REQUIRE(chip_no == max_idxs[(i >> 1)].second);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(max_score, tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(max_cfos, tolerance));
        }
        delete proc;
    }

    SECTION("Square Score") {
        constexpr float threshold = 220500.f;

        const vector<vector<float>> scores_out
            = load_data_matrix<float, float>(data_file, "score_sqr_raw_m10dB_w34_q64_N60_2pi_1_n30");

        vector<std::pair<int64_t, int64_t>> det_idxs(30, {0, 0});
        vector<std::pair<int64_t, int64_t>> max_idxs(30, {0, 0});
        for (unsigned u = 0; u < 30U; u++) {
            vector<int64_t> local_det_idx(p_omega);
            vector<float>   local_det_score(p_omega);
            vector<int64_t> local_max_idx(p_omega);
            vector<float>   local_max_score(p_omega);
            for (unsigned w = 0; w < p_omega; w++) {
                const float * begin = scores_out[w].data() + u * run_size;
                const float * end   = scores_out[w].data() + (u + 1) * run_size;

                const float * det_ptr = std::find_if(begin, end, [&](float value) { return value > threshold; });
                local_det_score[w]    = *det_ptr;
                local_det_idx[w]      = det_ptr - scores_out[w].data();

                const float * max_ptr = std::max_element(begin, end);
                local_max_score[w]    = *max_ptr;
                local_max_idx[w]      = max_ptr - scores_out[w].data();
            }

            // * Search the minimum detection index, and if there is several
            // * frequencies that exceed the threshold at the same time, the one
            // * corresponding to the highest score is kept.
            int64_t current_freq = 0;
            for (unsigned w = 0; w < p_omega; w++) {
                const int current_idx = local_det_idx[current_freq];
                const int local_idx   = local_det_idx[w];

                if (local_idx < current_idx) {
                    current_freq = w;
                } else if (local_idx == current_idx) {
                    current_freq = local_det_score[w] > local_det_score[current_freq] ? w : current_freq;

                    // * NOTE: Uncomment following lines to monitor this case (two frequency hypothesis exceed the threshold at the same time)
                    // printf("%3u | %3lu -- local %10d: %16.8f current %2d: %16.8f\n\n",
                    //        u, current_freq,
                    //        local_idx, local_det_score[w],
                    //        current_idx, local_det_score[current_freq]);
                }
            }
            const unsigned det_freq = current_freq;

            det_idxs[u] = {det_freq, local_det_idx[det_freq]};

            const float *  max_freq_score = std::max_element(local_max_score.data(), local_max_score.data() + p_omega);
            const unsigned max_freq       = max_freq_score - local_max_score.data();

            max_idxs[u] = {max_freq, local_max_idx[max_freq]};
        }

        detector_t * proc = new detector_t(
            pn.data(),
            threshold,
            step_denominator);

        for (unsigned i = 0; i < p_omega; i++) {
            const float proc_freq = proc->frequency_error(i) * 2 * q;
            const float ref_freq  = float(i) / step_denominator - float(p_omega - 1) / float(2 * step_denominator);
            REQUIRE_THAT(proc_freq, Catch::Matchers::WithinRel(ref_freq, 1e-6f));
        }

        REQUIRE(proc->threshold() == threshold);

        // * NOTE: Uncomment the following line and the two next groups to log the score
        // vector<vector<float>> result_scores(scores_out[0].size(), vector<float>(p_omega, 0));

        for (int64_t i = 0; i < int64_t(re_in.size()); i++) {
            proc->process_sqr(re_in[i], im_in[i], &running_state);

            // * NOTE: Uncomment the following line, the next group and the previous group to log the score
            // memcpy(result_scores[i].data(), running_state.scores, p_omega * sizeof(float));

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
            REQUIRE((running_state.max_score == running_state.scores[running_state.frequency_index] or (running_state.frame_detected)));
        }

        // * NOTE: Uncomment the following lines and the two previous groups to log the score
        // FILE * f = fopen("scores.dat", "w");
        // for (int j = 0; j < (int) result_scores.size(); j++) {
        //     fprintf(f, "%16.8f %16.8f %16.8f %16.8f\n", result_scores[j][0], result_scores[j][1], result_scores[j][2], result_scores[j][3]);
        // }
        // fclose(f);

        REQUIRE(results.size() == 60LU);

        constexpr float tolerance = float(5e-5);

        for (int64_t i = 0; i < int64_t(results.size()); i += 2) {
            const auto & result_i   = results[i];
            const auto & result_ip1 = results[i + 1];

            const int64_t chip_no = result_ip1.second - detector_t::window_size + result_ip1.first.chip_from_max;

            const float det_cfos = proc->frequency_error(det_idxs[i >> 1].first);
            const float max_cfos = proc->frequency_error(max_idxs[i >> 1].first);

            // * NOTE: If you want to explore the results, uncomment the 6 following lines
            // printf("-- Detection no %ld score exceeded the threshold of %5.1f at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e (%3d / %u) Hz.\n",
            //        i >> 1, proc->threshold(), result_i.second, result_i.first.max_score, result_i.first.frequency_offset,
            //        result_i.first.frequency_index, p_omega);
            // printf("-- Detection no %ld max has been found at chip no %ld,\n"
            //        "-- for a score of %16.8f and a frequency offset of %16.8e (%3d / %u) Hz.\n",
            //        i >> 1, chip_no, result_ip1.first.max_score, result_ip1.first.frequency_offset,
            //        result_i.first.frequency_index, p_omega);

            const float det_score = scores_out[det_idxs[i >> 1].first][det_idxs[i >> 1].second];
            const float max_score = scores_out[max_idxs[i >> 1].first][max_idxs[i >> 1].second];

            REQUIRE(result_i.first.frequency_index == det_idxs[i >> 1].first);
            REQUIRE(result_i.second == det_idxs[i >> 1].second);
            REQUIRE_THAT(result_i.first.max_score, Catch::Matchers::WithinRel(det_score, tolerance));
            REQUIRE_THAT(result_i.first.frequency_offset, Catch::Matchers::WithinRel(det_cfos, tolerance));
            REQUIRE(result_ip1.first.frequency_index == max_idxs[i >> 1].first);
            REQUIRE(chip_no == max_idxs[(i >> 1)].second);
            REQUIRE_THAT(result_ip1.first.max_score, Catch::Matchers::WithinRel(max_score, tolerance));
            REQUIRE_THAT(result_ip1.first.frequency_offset, Catch::Matchers::WithinRel(max_cfos, tolerance));
        }
        delete proc;
    }
    Mat_Close(data_file);
}
