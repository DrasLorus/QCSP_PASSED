#ifndef _QCSP_PASSED_CORR_ABS_MAX_GENERIC_HPP_
#define _QCSP_PASSED_CORR_ABS_MAX_GENERIC_HPP_ 1

#include <cstdint>
#include <vector>

namespace QCSP {
namespace StandaloneDetector {

class CCorrAbsMaxGeneric {
public:
    const unsigned q;    // Tq;
    const unsigned mask; // q - 1;

private:
    /**
     * @brief reference PN sequence
     *
     */
    const std::vector<float> pn; // [q];

    /**
     * @brief used PN sequence
     *
     */
    std::vector<float> rotating_pn; // [q];

    std::vector<float> re_corr_registers;  // [q];
    std::vector<float> im_corr_registers;  // [q];
    std::vector<float> abs_corr_registers; // [q];

public:
    const std::vector<float> & get_pn() const { return pn; }

    virtual float process(float re_in, float im_in);

    CCorrAbsMaxGeneric(const std::vector<float> & pn);

    virtual ~CCorrAbsMaxGeneric() = default;
};

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_CORR_ABS_MAX_GENERIC_HPP_
