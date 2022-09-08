#ifndef _QCSP_PASSED_MISC_HPP_
#define _QCSP_PASSED_MISC_HPP_ 1

#define _USE_MATH_DEFINES
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>

namespace QCSP {
namespace StandaloneDetector {

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

constexpr double pi      = (double) M_PI;
constexpr double two_pi  = 2. * pi;
constexpr double half_pi = pi / 2.;

constexpr float pi_f      = (float) M_PI;
constexpr float two_pi_f  = (float) (2. * pi);
constexpr float half_pi_f = (float) (pi / 2.);

template <typename T>
constexpr bool is_pow2(T value) {
    return value > 2
             ? is_pow2(value / 2)
             : (value == 2 ? true : 0);
}

template <uint64_t value>
constexpr uint8_t pow2_log2();

template <>
constexpr uint8_t pow2_log2<1>() {
    return 0;
}

template <>
constexpr uint8_t pow2_log2<2>() {
    return 1;
}

template <uint64_t value>
constexpr uint8_t pow2_log2() {
    static_assert(value > 0, "Value must be a power of 2.");
    static_assert(is_pow2(value), "Value must be a power of 2.");

    return value > 2
             ? pow2_log2<(value >> 1)>() + 1
             : 1U;
}

template <unsigned Tsize>
struct max_pow2;

template <>
struct max_pow2<0> {
    template <typename T>
    static constexpr inline T max(T *) {
        return 0.;
    }
};

template <>
struct max_pow2<1> {
    template <typename T>
    static constexpr inline T max(T * a) {
        return *a;
    }
};

template <>
struct max_pow2<2> {
    template <typename T>
    static constexpr inline T max(T * a) {
        return std::max(a[0], a[1]);
    }
};

template <unsigned Tsize>
struct max_pow2 {
    template <typename T>
    static constexpr inline T max(T * a) {
        static_assert(is_pow2(Tsize), "Size must be a power of 2.");
        // const T head = max_pow2<Tsize / 2>(a);
        // const T tail = max_pow2<Tsize / 2>(a + Tsize / 2);
        return std::max(max_pow2<Tsize / 2>::max(a), max_pow2<Tsize / 2>::max(a + Tsize / 2));
    }
};

template <typename T, T value>
constexpr bool is_even() {
    return value / 2U * 2U == value;
}

template <typename T, T value>
constexpr bool is_odd() {
    return not is_even<T, value>();
}

template <typename T, T saturation>
constexpr T low_sat(T a) {
    return std::max(a, saturation);
}

template <typename T>
constexpr T low_sat_1(T a) {
    return std::max(a, T(1));
}

template <typename T, T saturation>
constexpr T high_sat(T a) {
    return std::min(a, saturation);
}

constexpr float if_nan_0(float a) {
    return (std::isnan(a) ? 0 : a);
}

constexpr double if_nan_0(double a) {
    return (std::isnan(a) ? 0 : a);
}

} // namespace StandaloneDetector
} // namespace QCSP

#endif // _QCSP_PASSED_MISC_HPP_
