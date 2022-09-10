#ifndef _QCSP_PASSED_DETECTOR_HPP_
#define _QCSP_PASSED_DETECTOR_HPP_ 1

#include "./CScoreProcessor/CScoreProcessor.hpp"

#include <stdexcept>
#include <string>
#include <vector>

namespace QCSP {
namespace StandaloneDetector {

template <typename TScore, typename TFreq, unsigned Tp_omega>
struct DetectionState {
    static_assert(Tp_omega > 0U, "Can't have no hypothesis.");
    bool frame_detected;
    bool max_found;

    TScore scores[Tp_omega];
    TScore max_score;

    uint64_t chip_since_last_det;
    uint64_t chip_from_max;
    TFreq    frequency_offset; // The frequency offset corresponding to the maximum score
    uint32_t frequency_index;  // The index of the frequency offset in the frequency array
};

template <unsigned TFrameSize, unsigned Tq, unsigned Tp_omega, typename TIn_Type = float, bool normed = true, CorrelationEngineType variant = TIME_SLIDING>
class CDetector {
public:
    using state_t = DetectionState<TIn_Type, TIn_Type, Tp_omega>;

    virtual void process(TIn_Type re_in, TIn_Type im_in, state_t * state)     = 0;
    virtual void process_sqr(TIn_Type re_in, TIn_Type im_in, state_t * state) = 0;

    virtual ~CDetector() = default;
};

template <unsigned TFrameSize, unsigned Tq, unsigned Tp_omega, typename TIn_Type = float, bool normed = true, CorrelationEngineType variant = TIME_SLIDING>
class CDetectorSerial : public CDetector<TFrameSize, Tq, Tp_omega, TIn_Type, normed, variant> {
    static_assert(is_pow2(Tq), "q must be a power of 2");

public:
    static constexpr unsigned q       = Tq;
    static constexpr unsigned N       = TFrameSize;
    static constexpr unsigned p_omega = Tp_omega;

    using base = CDetector<N, q, p_omega, TIn_Type, normed>;
    using typename base::state_t;

    static constexpr uint64_t window_size = N * q;

private:
    using score_proc_t = CScoreProcessor<N, q, TIn_Type, normed, variant>;

    std::vector<score_proc_t> score_processors;

    TIn_Type frequency_errors[p_omega];

    const unsigned num_step;
    const unsigned den_step;

    const TIn_Type _symbol_rotation;
    const TIn_Type _rotation_step;

    const size_t          rotation_size;
    std::vector<TIn_Type> rotation_vect;
    size_t                rotation_increments[p_omega];
    size_t                rotation_counters[p_omega];

    TIn_Type _threshold;

    void update_state(const TIn_Type * __restrict scores, state_t * __restrict state) {
        const bool     frame_already_detected = state->frame_detected;
        const TIn_Type current_score          = state->max_score;
        const uint64_t current_idx            = state->chip_from_max;
        const TIn_Type current_cfos           = state->frequency_offset;
        const uint32_t current_cfid           = state->frequency_index;
        const uint64_t current_count          = state->chip_since_last_det + 1;

        const TIn_Type * max_score = std::max_element(scores, scores + p_omega);

        const bool higher_score   = current_score < *max_score;
        const bool over_threshold = *max_score > _threshold;
        const bool max_found      = current_count >= window_size;
        const bool relaxed        = current_count >= (window_size * 2) or current_count == 0;

        const uint32_t local_cfid = uint32_t(max_score - scores);
        const TIn_Type local_cfos = frequency_errors[local_cfid];

        const bool fetch_new_values = not(frame_already_detected and not(relaxed)) or (higher_score and not(max_found) and not(relaxed));

        const TIn_Type new_max_score = TIn_Type(fetch_new_values) * *max_score + TIn_Type(not fetch_new_values) * current_score;
        const TIn_Type new_cfos      = TIn_Type(fetch_new_values) * local_cfos + TIn_Type(not fetch_new_values) * current_cfos;
        const TIn_Type new_cfid      = TIn_Type(fetch_new_values) * local_cfid + TIn_Type(not fetch_new_values) * current_cfid;
        const uint64_t new_count     = current_count % (window_size * 2);

        memcpy(state->scores, scores, sizeof(TIn_Type) * p_omega);
        state->max_score           = new_max_score;
        state->frequency_offset    = new_cfos;
        state->frequency_index     = uint32_t(new_cfid);
        state->frame_detected      = (frame_already_detected and not(relaxed)) or over_threshold;
        state->max_found           = max_found and not(relaxed) and frame_already_detected;
        state->chip_since_last_det = uint64_t(frame_already_detected) * new_count;
        state->chip_from_max       = uint64_t(fetch_new_values) * new_count + uint64_t(not fetch_new_values) * current_idx;
    }

public:
    const TIn_Type * pn() const { return score_processors[0].get_pn(); }

    TIn_Type threshold() const noexcept { return _threshold; }
    TIn_Type symbol_rotation() const noexcept { return _symbol_rotation; }
    TIn_Type rotation_step() const noexcept { return _rotation_step; }

    const std::vector<TIn_Type> & primary_rotation() const noexcept { return rotation_vect; }

    TIn_Type frequency_error(unsigned n) const {
        if (n >= p_omega) {
            throw std::out_of_range("n must be below p_omega (= " + std::to_string(p_omega) + ")");
        }
        return frequency_errors[n];
    }

    virtual void process(TIn_Type re_in, TIn_Type im_in, state_t * state) override {
        TIn_Type scores[p_omega];

        for (unsigned u = 0; u < p_omega; u++) {
            const size_t current_counter = rotation_counters[u];
            const size_t local_counter   = current_counter << 1;

            const TIn_Type re_rotation = rotation_vect[local_counter];
            const TIn_Type im_rotation = rotation_vect[local_counter + 1];

            const TIn_Type local_re_in = re_in * re_rotation - im_in * im_rotation;
            const TIn_Type local_im_in = re_in * im_rotation + im_in * re_rotation;

            const TIn_Type score_u = score_processors[u].process(local_re_in, local_im_in);

            scores[u]            = score_u;
            rotation_counters[u] = (current_counter + rotation_increments[u]) % rotation_size;
        }

        update_state(scores, state);
    }

    virtual void process_sqr(TIn_Type re_in, TIn_Type im_in, state_t * state) override {
        TIn_Type scores[p_omega];

        for (unsigned u = 0; u < p_omega; u++) {
            const size_t current_counter = rotation_counters[u];
            const size_t local_counter   = current_counter << 1;

            const TIn_Type re_rotation = rotation_vect[local_counter];
            const TIn_Type im_rotation = rotation_vect[local_counter + 1];

            const TIn_Type local_re_in = re_in * re_rotation - im_in * im_rotation;
            const TIn_Type local_im_in = re_in * im_rotation + im_in * re_rotation;

            const TIn_Type score_u = score_processors[u].process_sqr(local_re_in, local_im_in);

            scores[u]            = score_u;
            rotation_counters[u] = (current_counter + rotation_increments[u]) % rotation_size;
        }

        update_state(scores, state);
    }

    template <typename Tpn>
    CDetectorSerial(Tpn * _pn, TIn_Type threshold, unsigned step_denominator, unsigned step_numerator = 1)
        : score_processors(p_omega, score_proc_t(_pn)),
          num_step(step_numerator * unsigned(p_omega > 1) * unsigned(step_denominator > 0)),
          den_step(step_denominator * unsigned(p_omega > 1)),
          _symbol_rotation(TIn_Type(pi_f * float(num_step) / low_sat_1(float(den_step)))),
          _rotation_step(_symbol_rotation / TIn_Type(q)),
          rotation_size(low_sat_1(2 * q * den_step)),
          rotation_vect(rotation_size * 2, 0),
          _threshold(threshold) {

        memset(rotation_counters, 0, p_omega * sizeof(size_t));

        std::vector<int> spanf(p_omega);
        // spanf = @(p_omega) -(p_omega - 1) : 2 : (p_omega - 1)
        const int span_root = -int(p_omega - 1);
        for (unsigned u = 0; u < p_omega; u++) {
            spanf[u] = span_root + int(u << 1);
        }

        for (unsigned u = 0; u < p_omega; u++) {
            frequency_errors[u] = rotation_step() / two_pi_f * TIn_Type(spanf[u]);
        }

        // rotation_vect = @(p_omega, den_step, q) exp(1i * pi / (q * den_step) .* (0 : (rotations_size - 1)))
        const double common_value = rotation_step();
        for (int i = 0; i < int(rotation_size); i++) {
            const int idx          = i << 1;
            rotation_vect[idx]     = (TIn_Type) cos(double(i) * common_value);
            rotation_vect[idx + 1] = (TIn_Type) sin(double(i) * common_value);
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

    virtual ~CDetectorSerial() = default;
};

template <unsigned TFrameSize, unsigned Tq, unsigned Tp_omega, typename TIn_Type = float, bool normed = true, CorrelationEngineType variant = TIME_SLIDING>
class CDetectorSerialMult : public CDetector<TFrameSize, Tq, Tp_omega, TIn_Type, normed, variant> {
    static_assert(is_pow2(Tq), "q must be a power of 2");

public:
    static constexpr unsigned q       = Tq;
    static constexpr unsigned N       = TFrameSize;
    static constexpr unsigned p_omega = Tp_omega;

    using base = CDetector<N, q, p_omega, TIn_Type, normed>;
    using typename base::state_t;

    static constexpr uint64_t window_size = N * q;

private:
    using score_proc_t = CScoreProcessor<N, q, TIn_Type, normed, variant>;

    std::vector<score_proc_t> score_processors;

    TIn_Type frequency_errors[p_omega];

    const unsigned num_step;
    const unsigned den_step;

    TIn_Type              _symbol_rotation;
    TIn_Type              _rotation_step;
    std::vector<TIn_Type> _root_rotations;
    std::vector<TIn_Type> _symbol_root_rotations;
    const size_t          rotation_size;
    size_t                rotation_counter;

    TIn_Type _threshold;

    void update_state(const TIn_Type * __restrict scores, state_t * __restrict state) {
        const bool     frame_already_detected = state->frame_detected;
        const TIn_Type current_score          = state->max_score;
        const uint64_t current_idx            = state->chip_from_max;
        const TIn_Type current_cfos           = state->frequency_offset;
        const uint32_t current_cfid           = state->frequency_index;
        const uint64_t current_count          = state->chip_since_last_det + 1;

        const TIn_Type * max_score = std::max_element(scores, scores + p_omega);

        const bool higher_score   = current_score < *max_score;
        const bool over_threshold = *max_score > _threshold;
        const bool max_found      = current_count >= window_size;
        const bool relaxed        = current_count >= (window_size * 2) or current_count == 0;

        const uint32_t local_cfid = uint32_t(max_score - scores);
        const TIn_Type local_cfos = frequency_errors[local_cfid];

        const bool fetch_new_values = not(frame_already_detected and not(relaxed)) or (higher_score and not(max_found) and not(relaxed));

        const TIn_Type new_max_score = TIn_Type(fetch_new_values) * *max_score + TIn_Type(not fetch_new_values) * current_score;
        const TIn_Type new_cfos      = TIn_Type(fetch_new_values) * local_cfos + TIn_Type(not fetch_new_values) * current_cfos;
        const TIn_Type new_cfid      = TIn_Type(fetch_new_values) * local_cfid + TIn_Type(not fetch_new_values) * current_cfid;
        const uint64_t new_count     = current_count % (window_size * 2);

        memcpy(state->scores, scores, sizeof(TIn_Type) * p_omega);
        state->max_score           = new_max_score;
        state->frequency_offset    = new_cfos;
        state->frequency_index     = uint32_t(new_cfid);
        state->frame_detected      = (frame_already_detected and not(relaxed)) or over_threshold;
        state->max_found           = max_found and not(relaxed) and frame_already_detected;
        state->chip_since_last_det = uint64_t(frame_already_detected) * new_count;
        state->chip_from_max       = uint64_t(fetch_new_values) * new_count + uint64_t(not fetch_new_values) * current_idx;
    }

public:
    const TIn_Type * pn() const { return score_processors[0].get_pn(); }

    TIn_Type threshold() const noexcept { return _threshold; }
    TIn_Type symbol_rotation() const noexcept { return _symbol_rotation; }
    TIn_Type rotation_step() const noexcept { return _rotation_step; }

    const std::vector<TIn_Type> & root_rotations() const noexcept { return _root_rotations; }
    const std::vector<TIn_Type> & symbol_root_rotations() const noexcept { return _symbol_root_rotations; }

    TIn_Type frequency_error(unsigned n) const {
        if (n >= p_omega) {
            throw std::out_of_range("n must be below p_omega (= " + std::to_string(p_omega) + ")");
        }
        return frequency_errors[n];
    }

    virtual void process(TIn_Type re_in, TIn_Type im_in, state_t * state) override {
        TIn_Type scores[p_omega];

        const size_t current_counter = rotation_counter;
        for (unsigned u = 0; u < p_omega; u++) {
            const TIn_Type common_value = _root_rotations[u] * TIn_Type(current_counter);

            const TIn_Type re_rotation = (TIn_Type) std::cos(common_value);
            const TIn_Type im_rotation = (TIn_Type) std::sin(common_value);

            const TIn_Type local_re_in = re_in * re_rotation - im_in * im_rotation;
            const TIn_Type local_im_in = re_in * im_rotation + im_in * re_rotation;
            const TIn_Type score_u     = score_processors[u].process(local_re_in, local_im_in);

            scores[u] = score_u;
        }
        rotation_counter = (current_counter + 1) % rotation_size;

        update_state(scores, state);
    }

    virtual void process_sqr(TIn_Type re_in, TIn_Type im_in, state_t * state) override {
        TIn_Type scores[p_omega];

        const size_t current_counter = rotation_counter;
        for (unsigned u = 0; u < p_omega; u++) {
            const TIn_Type common_value = _root_rotations[u] * TIn_Type(current_counter);

            const TIn_Type re_rotation = (TIn_Type) std::cos(common_value);
            const TIn_Type im_rotation = (TIn_Type) std::sin(common_value);

            const TIn_Type local_re_in = re_in * re_rotation - im_in * im_rotation;
            const TIn_Type local_im_in = re_in * im_rotation + im_in * re_rotation;

            const TIn_Type score_u = score_processors[u].process_sqr(local_re_in, local_im_in);

            scores[u] = score_u;
        }
        rotation_counter = (current_counter + 1) % rotation_size;

        update_state(scores, state);
    }

    template <typename Tpn>
    CDetectorSerialMult(Tpn * _pn, TIn_Type threshold, unsigned step_denominator, unsigned step_numerator = 1)
        : score_processors(p_omega, score_proc_t(_pn)),
          num_step(step_numerator * unsigned(p_omega > 1) * unsigned(step_denominator > 0)),
          den_step(step_denominator * unsigned(p_omega > 1)),
          _symbol_rotation(TIn_Type(pi_f * float(num_step) / low_sat_1(float(den_step)))),
          _rotation_step(_symbol_rotation / TIn_Type(q)),
          _root_rotations(p_omega, TIn_Type(0)),
          _symbol_root_rotations(p_omega, TIn_Type(0)),
          rotation_size(low_sat_1(2 * q * den_step)),
          rotation_counter(0),
          _threshold(threshold) {

        const TIn_Type   common_value = rotation_step();
        std::vector<int> spanf(p_omega);
        // spanf = @(p_omega) -(p_omega - 1) : 2 : (p_omega - 1)
        const int span_root = -int(p_omega - 1);
        for (unsigned u = 0; u < p_omega; u++) {
            spanf[u]                  = span_root + int(u << 1);
            _root_rotations[u]        = TIn_Type(spanf[u]) * common_value;
            _symbol_root_rotations[u] = _root_rotations[u] * TIn_Type(q);
        }

        for (unsigned u = 0; u < p_omega; u++) {
            frequency_errors[u] = rotation_step() / two_pi_f * TIn_Type(spanf[u]);
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

    virtual ~CDetectorSerialMult() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_DETECTOR_HPP_
