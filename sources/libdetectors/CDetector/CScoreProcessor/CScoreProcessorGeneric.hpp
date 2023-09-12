#ifndef _QCSP_PASSED_SCORE_PROCESSOR_GENERIC_HPP_
#define _QCSP_PASSED_SCORE_PROCESSOR_GENERIC_HPP_ 1

#include <cstdint>
#include <vector>

#include "./CCorrelationEngine/CCorrelationEngineGeneric.hpp"
#include "./CNorm/CNormGeneric.hpp"
#include "./CScoreAccumulator/CScoreAccumulatorGeneric.hpp"

namespace QCSP {
namespace StandaloneDetector {

class CScoreProcessorGeneric {
public:
    const unsigned q;
    const unsigned N;
    const bool     normed;

private:
    CCorrelationEngineGeneric corr_engine;

    CNormGeneric             norm_proc;
    CScoreAccumulatorGeneric score_accumulator;

    /**
     * @brief Conditional pointer to allow runtime normalization
     *
     * @details conditional_norm points allow to either notmalize, or not.
     */
    float (*conditional_norm)(float, float);

public:
    const std::vector<float> & get_pn() const;

    virtual float process(float re_in, float im_in);

    virtual float process_sqr(float re_in, float im_in);

    CScoreProcessorGeneric(const std::vector<float> & _pn, unsigned _N, bool _normed = true);

    virtual ~CScoreProcessorGeneric() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_SCORE_PROCESSOR_GENERIC_HPP_
