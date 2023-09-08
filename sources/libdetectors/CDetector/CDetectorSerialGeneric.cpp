#include <algorithm>
#include <cstdint>
#include <cstring>
#include <stdexcept>
#include <string>
#include <vector>

#include "Miscellanous/misc.hpp"

#include "./CDetectionStateGeneric.hpp"
#include "./CDetectorGenericInterface.hpp"
#include "./CDetectorSerialGeneric.hpp"

void QCSP::StandaloneDetector::CDetectorSerialGeneric::update_state(const float * __restrict scores, state_t * __restrict state) const {
    const bool     frame_already_detected = state->frame_detected;
    const float    current_score          = state->max_score;
    const uint64_t current_idx            = state->chip_from_max;
    const float    current_cfos           = state->frequency_offset;
    const uint32_t current_cfid           = state->frequency_index;
    const uint64_t current_count          = state->chip_since_last_det + 1;

    const float * max_score = std::max_element(scores, scores + p_omega);

    const bool higher_score   = current_score < *max_score;
    const bool over_threshold = *max_score > _threshold;
    const bool max_found      = current_count >= window_size;
    const bool relaxed        = current_count >= (window_size * 2) or current_count == 0;

    const uint32_t local_cfid = uint32_t(max_score - scores);
    const float    local_cfos = frequency_errors[local_cfid];

    const bool fetch_new_values = not(frame_already_detected and not(relaxed)) or (higher_score and not(max_found) and not(relaxed));

    const float    new_max_score = float(fetch_new_values) * *max_score + float(not fetch_new_values) * current_score;
    const float    new_cfos      = float(fetch_new_values) * local_cfos + float(not fetch_new_values) * current_cfos;
    const float    new_cfid      = float(fetch_new_values) * local_cfid + float(not fetch_new_values) * current_cfid;
    const uint64_t new_count     = current_count % (window_size * 2);

    std::memcpy(state->scores.data(), scores, sizeof(float) * p_omega);
    state->max_score           = new_max_score;
    state->frequency_offset    = new_cfos;
    state->frequency_index     = uint32_t(new_cfid);
    state->frame_detected      = (frame_already_detected and not(relaxed)) or over_threshold;
    state->max_found           = max_found and not(relaxed) and frame_already_detected;
    state->chip_since_last_det = uint64_t(frame_already_detected) * new_count;
    state->chip_from_max       = uint64_t(fetch_new_values) * new_count + uint64_t(not fetch_new_values) * current_idx;
}

void QCSP::StandaloneDetector::CDetectorSerialGeneric::process(float re_in, float im_in, state_t * state) {
    float * scores = local_scores.data();

    for (unsigned u = 0; u < p_omega; u++) {
        const size_t current_counter = rotation_counters[u];
        const size_t local_counter   = current_counter << 1;

        const float re_rotation = rotation_vect[local_counter];
        const float im_rotation = rotation_vect[local_counter + 1];

        const float local_re_in = re_in * re_rotation - im_in * im_rotation;
        const float local_im_in = re_in * im_rotation + im_in * re_rotation;

        const float score_u = score_processors[u].process(local_re_in, local_im_in);

        scores[u]            = score_u;
        rotation_counters[u] = (current_counter + rotation_increments[u]) % rotation_size;
    }

    update_state(scores, state);
}

void QCSP::StandaloneDetector::CDetectorSerialGeneric::process_sqr(float re_in, float im_in, state_t * state) {
    float * scores = local_scores.data();

    for (unsigned u = 0; u < p_omega; u++) {
        const size_t current_counter = rotation_counters[u];
        const size_t local_counter   = current_counter << 1;

        const float re_rotation = rotation_vect[local_counter];
        const float im_rotation = rotation_vect[local_counter + 1];

        const float local_re_in = re_in * re_rotation - im_in * im_rotation;
        const float local_im_in = re_in * im_rotation + im_in * re_rotation;

        const float score_u = score_processors[u].process_sqr(local_re_in, local_im_in); //! FIXME

        scores[u]            = score_u;
        rotation_counters[u] = (current_counter + rotation_increments[u]) % rotation_size;
    }

    update_state(scores, state);
}

QCSP::StandaloneDetector::CDetectorSerialGeneric::CDetectorSerialGeneric(
    const std::vector<float> & _pn,
    uint32_t                   _N,
    uint32_t                   _p_omega,
    float                      _threshold,
    unsigned                   _step_denominator,
    unsigned                   _step_numerator,
    bool                       _normed)
    : q(_pn.size()),
      N(_N),
      p_omega(_p_omega),
      normed(_normed),
      window_size(N * q),
      score_processors(p_omega, score_proc_t(_pn, N, normed)), //! FIXME
      local_scores(p_omega, 0),
      frequency_errors(p_omega, 0),
      num_step(_step_numerator * unsigned(p_omega > 1) * unsigned(_step_denominator > 0)),
      den_step(_step_denominator * unsigned(p_omega > 1)),
      _symbol_rotation(float(pi_f * float(num_step) / low_sat_1(float(den_step)))),
      _rotation_step(_symbol_rotation / float(q)),
      rotation_size(low_sat_1(2 * q * den_step)),
      rotation_vect(rotation_size * 2, 0),
      rotation_increments(p_omega, 0),
      rotation_counters(p_omega, 0),
      _threshold(_threshold) {

    std::vector<int> spanf(p_omega);
    // spanf = @(p_omega) -(p_omega - 1) : 2 : (p_omega - 1)
    const int span_root = -int(p_omega - 1);
    for (unsigned u = 0; u < p_omega; u++) {
        spanf[u] = span_root + int(u << 1);
    }

    for (unsigned u = 0; u < p_omega; u++) {
        frequency_errors[u] = rotation_step() / two_pi_f * float(spanf[u]);
    }

    // rotation_vect = @(p_omega, den_step, q) exp(1i * pi / (q * den_step) .* (0 : (rotations_size - 1)))
    const double common_value = rotation_step();
    for (int i = 0; i < int(rotation_size); i++) {
        const int idx          = i << 1;
        rotation_vect[idx]     = (float) cos(double(i) * common_value);
        rotation_vect[idx + 1] = (float) sin(double(i) * common_value);
    }

    for (int u = 0; u < int(p_omega >> 1); u++) {
        rotation_increments[u] = size_t(spanf[u] + int(rotation_size));
    }
    for (unsigned u = (p_omega >> 1); u < p_omega; u++) {
        rotation_increments[u] = spanf[u];
    }

    // * NOTE: Uncomment following lines to monitor spanf and rotation_increments
    // for (unsigned u = 0; u < p_omega; u++) {
    //     printf("%3d ", spanf[u]);
    // }
    // printf("\n");
    // for (unsigned u = 0; u < p_omega; u++) {
    //     printf("%3lu ", rotation_increments[u]);
    // }
    // printf("\n\n");
}
