#include <cmath>
#include <cstdint>
#include <cstring>
#include <memory>
#include <vector>

#include "./CScoreProcessorGeneric.hpp"

namespace {

float normalize(float x, float y) { return x / y; }
float identity(float x, float) { return x; }

} // namespace

const std::vector<float> & QCSP::StandaloneDetector::CScoreProcessorGeneric::get_pn() const { return corr_engine.get_pn(); }

float QCSP::StandaloneDetector::CScoreProcessorGeneric::process(float re_in, float im_in) {
    const float norm_unsafe = norm_proc.process(re_in, im_in);
    const float norm_value  = norm_unsafe < float(1e-8)
                                ? float(1)
                                : norm_unsafe;

    const float cabs_max = std::sqrt(conditional_norm(corr_engine.process(re_in, im_in), norm_value));

    return score_accumulator.process(cabs_max);
}

float QCSP::StandaloneDetector::CScoreProcessorGeneric::process_sqr(float re_in, float im_in) {
    const float norm_unsafe = norm_proc.process(re_in, im_in);
    const float norm_value  = norm_unsafe < float(1e-8)
                                ? float(1)
                                : norm_unsafe;

    const float cabs_max = conditional_norm(corr_engine.process(re_in, im_in), norm_value);

    return score_accumulator.process(cabs_max);
}

QCSP::StandaloneDetector::CScoreProcessorGeneric::CScoreProcessorGeneric(const std::vector<float> & _pn, unsigned _N, bool _normed)
    : q(_pn.size()),
      N(_N),
      normed(_normed),
      corr_engine(_pn),
      norm_proc(q),
      score_accumulator(q, N),
      conditional_norm(_normed ? &normalize : &identity) {
}
