#ifndef _QCSP_PASSED_CGENERATOR_HPP_
#define _QCSP_PASSED_CGENERATOR_HPP_

#include <complex>
#include <cstdio>
#include <random>
#include <string>
#include <vector>

#include "CGenerator/nbldpc_functions/tools.hpp"
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

    using complexf = std::complex<float>;

    struct frame_param_t {
        int32_t delay;
        float   rotation;
        float   phase;
    };

private:
    //* For NB-LDPC Generation
    code_t  code;
    table_t table;

    string _alist_file;

    //* For QCSP-specific tasks
    vector<bool> _pn;
    vector<bool> _ovmod;

    //* For random generation
    const float _snr;
    const float _sigma;
    const float _rotation_span;

    std::mt19937_64 gen_message;
    std::mt19937_64 gen_delay;
    std::mt19937_64 gen_phase;
    std::mt19937_64 gen_rotation;
    std::mt19937_64 gen_noise;

    std::uniform_int_distribution<int32_t> distr_message;
    std::uniform_int_distribution<int32_t> distr_delay;
    std::uniform_real_distribution<float>  distr_phase;
    std::uniform_real_distribution<float>  distr_rotation;
    std::normal_distribution<float>        distr_noise;

    //** Private methods
    void message_generation(vector<int32_t> & message);

    void nbldpc_encode(const vector<int32_t> & message, vector<int32_t> & codeword) const;
    void ccsk_modulate(const vector<int32_t> & codeword, vector<bool> & ccsk_frame) const;
    void overmodulate(const vector<bool> & ccsk_frame, vector<bool> & qcsp_frame) const;

    void add_delay(const vector<bool> &  qcsp_frame,
                   const frame_param_t & param,
                   vector<bool> &        delayed_frame) const; // Add delay to the frame, if any
    void add_bpsk_phase_rotation(const vector<bool> &  delayed_frame,
                                 const frame_param_t & param,
                                 vector<complexf> &    rotated_frame) const; // Add BPSK + add phase and rotation if any

    void add_noise(const vector<complexf> & rotated_sequence, vector<complexf> & noisy_sequence); // Generate and add white Gaussian noise

public:
    inline int32_t N() const { return code.N; }
    inline int32_t K() const { return code.K; }
    inline int32_t M() const { return code.M; }
    inline int32_t q() const { return code.q; }
    inline int32_t p() const { return code.p; }
    inline int32_t frame_size() const { return code.N * code.q; }
    inline int32_t run_length() const { return code.N * code.q * 5; }

    float snr() const { return _snr; };
    float sigma() const { return _sigma; };
    float rotation_span() const { return _rotation_span; };

    /**
     * @brief Construct a new CGenerator object
     *
     * @param alist_file Path to the alist file of the chosen code.
     * @param pn CCSK sequence. Its size must match the order of the code. The behavior is undefined otherwise.
     * @param ovmod Overmodulation sequence. Its size must match the size of the codeword. The behavior is undefined otherwise.
     * @param snr Signal to Noise Ratio targeted.
     * @param rotation_span The range of rotation. The resulting rotation will be between \f$\frac{-rotation_span}{2}\f$ and \f$\frac{rotation_span}{2}.
     */
    CGenerator(const string & alist_file, const vector<float> & pn, const vector<float> & ovmod, float snr, float rotation_span)
        : _alist_file(alist_file),
          _snr(snr),
          _sigma(sqrtf(powf(10.f, (-snr / 10.f)) / 2.f)),
          _rotation_span(rotation_span),
          gen_message(std::random_device()()),
          gen_delay(std::random_device()()),
          gen_phase(std::random_device()()),
          gen_rotation(std::random_device()()),
          gen_noise(std::random_device()()),
          distr_message(),
          distr_delay(),
          distr_phase(-pi_f, pi_f),
          distr_rotation(-rotation_span / 2.f, rotation_span / 2.f),
          distr_noise(0, _sigma) {

        //* First, load the alist !!
        LoadCode(_alist_file.c_str(), code);
        LoadTables(table, code.q, code.p);
        GaussianElimination(table, code);

        _pn = vector<bool>(code.q);
        for (int i = 0; i < code.q; i++) {
            _pn[i] = pn[i] > 0;
        }

        _ovmod = vector<bool>(code.N);
        for (int i = 0; i < code.N; i++) {
            _ovmod[i] = ovmod[i] > 0;
        }

        distr_message.param(std::uniform_int_distribution<int32_t>::param_type {0, code.q - 1});
        distr_delay.param(std::uniform_int_distribution<int32_t>::param_type {0, code.q - 1});
    }

    template <bool noise_only = false>
    frame_param_t generate(vector<complexf> & noisy_sequence);

    virtual ~CGenerator() = default;
};

} // namespace Simulators

} // namespace StandaloneDetector

} // namespace QCSP

#endif // _QCSP_PASSED_CGENERATOR_HPP_
