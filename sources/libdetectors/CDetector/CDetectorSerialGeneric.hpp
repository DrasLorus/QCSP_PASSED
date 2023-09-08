#ifndef _QCSP_PASSED_DETECTOR_SERIAL_GENERIC_HPP_
#define _QCSP_PASSED_DETECTOR_SERIAL_GENERIC_HPP_ 1

#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

#include "./CDetectionStateGeneric.hpp"
#include "./CDetectorGenericInterface.hpp"

#include "./CScoreProcessor/CScoreProcessorGeneric.hpp"

namespace QCSP {
namespace StandaloneDetector {

class CDetectorSerialGeneric : CDetectorGeneric {
public:
    const unsigned q;
    const unsigned N;
    const unsigned p_omega;
    const bool     normed;

    using state_t = DetectionStateGeneric;

    const uint64_t window_size; // = N * q;

private:
    using score_proc_t = CScoreProcessorGeneric;

    std::vector<score_proc_t> score_processors;
    std::vector<float>        local_scores;

    std::vector<float> frequency_errors; // [p_omega];

    const unsigned num_step;
    const unsigned den_step;

    const float _symbol_rotation;
    const float _rotation_step;

    const size_t        rotation_size;
    std::vector<float>  rotation_vect;
    std::vector<size_t> rotation_increments; // [p_omega];
    std::vector<size_t> rotation_counters;   // [p_omega];

    float _threshold;

    void update_state(const float * __restrict scores, state_t * __restrict state) const;

public:
    const std::vector<float> & pn() const { return score_processors[0].get_pn(); }

    float threshold() const noexcept { return _threshold; }
    float symbol_rotation() const noexcept { return _symbol_rotation; }
    float rotation_step() const noexcept { return _rotation_step; }

    const std::vector<float> & primary_rotation() const noexcept { return rotation_vect; }

    float frequency_error(unsigned n) const {
        if (n >= p_omega) {
            throw std::out_of_range("n must be below p_omega (= " + std::to_string(p_omega) + ")");
        }
        return frequency_errors[n];
    }

    virtual void process(float re_in, float im_in, state_t * state) override;

    virtual void process_sqr(float re_in, float im_in, state_t * state) override;

    CDetectorSerialGeneric(
        const std::vector<float> & pn,
        uint32_t                   N,
        uint32_t                   p_omega,
        float                      threshold,
        unsigned                   step_denominator = 1,
        unsigned                   step_numerator   = 1,
        bool                       normed           = true);

    virtual ~CDetectorSerialGeneric() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif
