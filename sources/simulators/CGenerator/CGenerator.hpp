#ifndef _QCSP_PASSED_CGENERATOR_HPP_
#define _QCSP_PASSED_CGENERATOR_HPP_

#include <string>
#include <vector>

#include "Miscellanous/misc.hpp"

#include "./nbldpc_functions/nbldpc_wrapper.hpp"

namespace QCSP {

namespace StandaloneDetector {

namespace Simulators {

class CGenerator {
public:
    using string = std::string;

    template <class T>
    using vector = std::vector<T>;

private:
    NbldpcWrapper::code_t  code;
    NbldpcWrapper::table_t table;

    vector<float> _pn;
    vector<float> _ovmod;

    string _alist_file;

    const float _snr;
    const float _sigma;
    const float _rotation_span;

    void message_generation(vector<int32_t> & message);
    void nbldpc_encode(const vector<int32_t> & message, vector<int32_t> & codeword);
    void ccsk_modulate(const vector<int32_t> & codeword, vector<int32_t> & ccsk_frame);
    void overmodulate(const vector<int32_t> & ccsk_frame, vector<int32_t> & qcsp_frame);

    void add_delay(const vector<int32_t> & noisy_seq, vector<int32_t> & delayed_frame);          // Add delay to the frame
    void add_phase_rotation(const vector<int32_t> & delayed_frame, vector<float> rotated_frame); // Add phase and rotation if any
    void add_noise(const vector<float> rotated_sequence, vector<float> noisy_sequence);          // Generate and add white Gaussian noise

public:
    template <class Tpn, class Tovmod>
    CGenerator(const string & alist_file, Tpn * pn, Tovmod * ovmod, float snr, float rotation_span)
        : _alist_file(alist_file),
          _snr(snr),
          _sigma(sqrtf(powf(10.f, (-snr / 10.f)) / 2.f)),
          _rotation_span(rotation_span) {

        NbldpcWrapper::LoadCode(_alist_file.c_str(), code);
        NbldpcWrapper::LoadTables(table, code.q, code.p);

        _pn    = vector<float>(pn, pn + code.q);
        _ovmod = vector<float>(ovmod, ovmod + code.N);
    }

    virtual ~CGenerator() = default;
};

} // namespace Simulators

} // namespace StandaloneDetector

} // namespace QCSP

#endif // _QCSP_PASSED_CGENERATOR_HPP_
