#ifndef _QCSP_METATYPES_HPP_
#define _QCSP_METATYPES_HPP_ 1

#include <cstdint>
#include <type_traits>

namespace QCSP {

namespace metatypes {

enum deduction_method : bool {
    FROM_MAX  = false,
    FROM_BITS = true
};

enum signedness : bool {
    UNSIGNED = false,
    SIGNED   = true
};

enum types : uint8_t {
    u8,
    u16,
    u32,
    u64,
    i8,
    i16,
    i32,
    i64
};

template <uint64_t Tvalue, enum signedness Tsigned, enum deduction_method method>
static constexpr enum types deduce_size() {
    if constexpr (method == FROM_MAX) {
        if constexpr (Tsigned == UNSIGNED) {
            if constexpr (Tvalue <= UINT8_MAX) {
                return u8;
            }
            if constexpr ((UINT8_MAX < Tvalue) and (Tvalue < UINT16_MAX)) {
                return u16;
            }
            if constexpr ((UINT16_MAX <= Tvalue) and (Tvalue < UINT32_MAX)) {
                return u32;
            }
            if constexpr ((UINT32_MAX <= Tvalue) && (Tvalue < UINT64_MAX)) {
                return u64;
            }
        }
        if constexpr (Tsigned == SIGNED) {
            if constexpr (Tvalue - INT8_MIN <= UINT8_MAX) {
                return i8;
            }
            if constexpr ((UINT8_MAX < int64_t(Tvalue)) and (int64_t(Tvalue) < UINT16_MAX)) {
                return i16;
            }
            if constexpr ((UINT16_MAX <= int64_t(Tvalue)) and (int64_t(Tvalue) < UINT32_MAX)) {
                return i32;
            }
            if constexpr ((UINT32_MAX <= int64_t(Tvalue)) && (int64_t(Tvalue) < UINT64_MAX)) {
                return i64;
            }
        }
    } else /*if constexpr (method != FROM_MAX) */ {
        if constexpr (Tsigned == UNSIGNED) {
            if constexpr (Tvalue <= 8) {
                return u8;
            }
            if constexpr (Tvalue <= 16) {
                return u16;
            }
            if constexpr (Tvalue <= 32) {
                return u32;
            }
            if constexpr (Tvalue <= 64) {
                return u64;
            }
        }
        if constexpr (Tsigned == SIGNED) {
            if constexpr (Tvalue <= 8) {
                return i8;
            }
            if constexpr (Tvalue <= 16) {
                return i16;
            }
            if constexpr (Tvalue <= 32) {
                return i32;
            }
            if constexpr (Tvalue <= 64) {
                return i64;
            }
        }
    }
}

template <enum types condition>
struct sintxx_gen {};

template <>
struct sintxx_gen<u8> {
    using type = uint8_t;
    template <class T>
    constexpr type operator()(T v) {
        return type(v);
    }
};

template <>
struct sintxx_gen<u16> {
    using type = uint16_t;
    template <class T>
    constexpr type operator()(T v) {
        return type(v);
    }
};

template <>
struct sintxx_gen<u32> {
    using type = uint32_t;
    template <class T>
    constexpr type operator()(T v) {
        return type(v);
    }
};

template <>
struct sintxx_gen<u64> {
    using type = uint64_t;
    template <class T>
    constexpr type operator()(T v) {
        return type(v);
    }
};

template <>
struct sintxx_gen<i8> {
    using type = int8_t;
    template <class T>
    constexpr type operator()(T v) {
        return type(v);
    }
};

template <>
struct sintxx_gen<i16> {
    using type = int16_t;
    template <class T>
    constexpr type operator()(T v) {
        return type(v);
    }
};

template <>
struct sintxx_gen<i32> {
    using type = int32_t;
    template <class T>
    constexpr type operator()(T v) {
        return type(v);
    }
};

template <>
struct sintxx_gen<i64> {
    using type = int64_t;
    template <class T>
    constexpr type operator()(T v) {
        return type(v);
    }
};

} // namespace metatypes

template <uint64_t Tvalue, metatypes::signedness Tsigned, enum metatypes::deduction_method method = metatypes::FROM_MAX>
using aintxx_t = typename metatypes::sintxx_gen<metatypes::deduce_size<Tvalue, Tsigned, method>()>::type;

template <uint64_t Tvalue, enum metatypes::deduction_method method = metatypes::FROM_MAX>
using uintxx_t = aintxx_t<Tvalue, metatypes::UNSIGNED, method>;

template <uint64_t Tvalue, enum metatypes::deduction_method method = metatypes::FROM_MAX>
using intxx_t = aintxx_t<Tvalue, metatypes::SIGNED, method>;

template <uint64_t Tvalue, metatypes::signedness Tsigned, enum metatypes::deduction_method method = metatypes::FROM_MAX>
constexpr aintxx_t<Tvalue, Tsigned, method> make_aintxx_t(aintxx_t<Tvalue, Tsigned, method> value = aintxx_t<Tvalue, Tsigned, method>(Tvalue)) {
    return aintxx_t<Tvalue, Tsigned, method>(value);
}

template <uint64_t Tvalue, enum metatypes::deduction_method method = metatypes::FROM_MAX>
constexpr uintxx_t<Tvalue, method> make_uintxx_t(uintxx_t<Tvalue, method> value = uintxx_t<Tvalue, method>(Tvalue)) {
    return make_aintxx_t<Tvalue, metatypes::UNSIGNED, method>(value);
}

template <uint64_t Tvalue, enum metatypes::deduction_method method = metatypes::FROM_MAX>
constexpr intxx_t<Tvalue, method> make_intxx_t(intxx_t<Tvalue, method> value = intxx_t<Tvalue, method>(Tvalue)) {
    return make_aintxx_t<Tvalue, metatypes::SIGNED, method>(value);
}

// template <uint64_t Tvalue, enum metatypes::deduction_method method = metatypes::FROM_MAX>
// using intxx_t = typename metatypes::sintxx_gen<metatypes::deduce_size<Tvalue, metatypes::SIGNED, method>()>::type;

// template <uint64_t Tvalue, metatypes::signedness Tsigned, enum metatypes::deduction_method method = metatypes::FROM_MAX>
// using aintxx_t = typename metatypes::sintxx_gen<metatypes::deduce_size<Tvalue, Tsigned, method>()>::type;

template <uint64_t in_bits, uint64_t truncated_bits, metatypes::signedness Tsigned = metatypes::UNSIGNED, bool lsb = true>
struct truncation {
    static_assert(in_bits >= truncated_bits, "Truncation cannot exceeds the bit-width.");
    static constexpr uint64_t out_bits = in_bits - truncated_bits;
    using in_t                         = aintxx_t<in_bits, Tsigned, metatypes::FROM_BITS>;
    using out_t                        = aintxx_t<out_bits, Tsigned, metatypes::FROM_BITS>;

    static constexpr out_t apply(in_t value) {
        if constexpr (lsb) {
            return out_t(value >> truncated_bits);
        } else {
            return out_t(value & ((1ULL << (out_bits)) - 1U));
        }
    }
};

template <uint64_t in_bits, uint64_t saturated_bits, metatypes::signedness Tsigned = metatypes::UNSIGNED>
struct saturation {
    static_assert(in_bits >= saturated_bits, "saturation cannot exceeds the bit-width.");
    static constexpr uint64_t out_bits = in_bits - saturated_bits;
    using in_t                         = aintxx_t<in_bits, Tsigned, metatypes::FROM_BITS>;
    using out_t                        = aintxx_t<out_bits, Tsigned, metatypes::FROM_BITS>;

    static constexpr out_t max_raw_value = out_t((1ULL << (out_bits - uint64_t(Tsigned))) - 1);
    static constexpr out_t min_raw_value = out_t(-(1ULL << (out_bits - uint64_t(Tsigned)))) * out_t(Tsigned);

    static constexpr out_t apply(in_t value) {
        if constexpr (Tsigned == metatypes::UNSIGNED) {
            return max_raw_value < value
                     ? max_raw_value
                     : out_t(value);
        } else {
            return max_raw_value < value
                     ? max_raw_value
                 : value < min_raw_value
                     ? min_raw_value
                     : out_t(value);
        }
    }
};


} // namespace QCSP

#endif // _QCSP_METATYPES_HPP_
