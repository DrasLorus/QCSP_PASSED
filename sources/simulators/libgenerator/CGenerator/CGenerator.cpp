#include "./CGenerator.hpp"
#include <algorithm>

#include <fstream>

using namespace QCSP::StandaloneDetector::Simulators;

void CGenerator::message_generation(vector<int32_t> & message) {
    for (int32_t i = 0; i < K(); i++) {
        message[i] = distr_message(gen_message);
    }
}

void CGenerator::nbldpc_encode(const vector<int32_t> & message, vector<int32_t> & codeword) const {
    Encoding(code, table, message, codeword);
}

void CGenerator::ccsk_modulate(const vector<int32_t> & codeword, vector<bool> & ccsk_frame) const {
    const auto iter_begin = _pn.cbegin();
    const auto iter_end   = _pn.cend();
    for (int32_t idx_s = 0; idx_s < N(); idx_s++) {
        const int32_t symbol      = codeword[idx_s];
        const auto    iter_middle = _pn.cbegin() + symbol;

        auto iter_dest = ccsk_frame.begin() + idx_s * q();

        // PN Rotation
        // do: [begin, ···, middle, ···, end] → [middle, ···, end, begin, ···, middle - 1]
        std::rotate_copy(iter_begin, iter_middle, iter_end, iter_dest);
    }
}

void CGenerator::overmodulate(const vector<bool> & ccsk_frame, vector<bool> & qcsp_frame) const {
    for (int32_t i = 0; i < N(); i++) {
        const bool over_symbol = _ovmod[i];

        const auto iter_begin = ccsk_frame.cbegin() + (i << p());
        const auto iter_end   = ccsk_frame.cbegin() + ((i + 1) << p());

        auto iter_dest = qcsp_frame.begin() + (i << p());

        std::transform(iter_begin, iter_end, iter_dest, [over_symbol](bool b) { return !(b ^ over_symbol); });
    }
}

void CGenerator::add_delay(const vector<bool> & qcsp_frame, const frame_param_t & param, vector<bool> & delayed_frame) const {
    const int32_t delay = param.delay;

    std::fill(delayed_frame.begin(), delayed_frame.end(), 0);

    const auto iter_dest = delayed_frame.begin() + frame_size() * 2 + delay;

    std::copy(qcsp_frame.cbegin(), qcsp_frame.cend(), iter_dest);
} // Add delay to the frame

void CGenerator::add_bpsk_phase_rotation(const vector<bool> & delayed_frame, const frame_param_t & param, vector<complexf> & rotated_frame) const {
    const float   phase    = param.phase;
    const float   rotation = param.rotation;
    const int32_t delay    = param.delay;

    for (int i = frame_size() * 2 + delay; i < frame_size() * 3 + delay; i++) {
        const double local_rotation = double(rotation) * double(i) / double(q()) + double(phase);

        rotated_frame[i] = complexf(std::exp(std::complex<double>(0., local_rotation)) * -double(int32_t(delayed_frame[i]) * 2 - 1));
    }
} // Add BPSK + add phase and rotation if any

void CGenerator::add_noise(const vector<complexf> & rotated_frame, vector<complexf> & noisy_sequence) {
    for (int32_t i = 0; i < run_length(); i++) {
        noisy_sequence[i] = rotated_frame[i] + complexf(distr_noise(gen_noise), distr_noise(gen_noise));
    }
}

template <>
CGenerator::frame_param_t CGenerator::generate<false>(vector<complexf> & noisy_sequence) {
    const frame_param_t run_parameters {distr_delay(gen_delay), distr_rotation(gen_rotation), distr_phase(gen_phase)};

    vector<int32_t> message(K()), codeword(N());
    vector<bool>    ccsk_frame(frame_size()), qcsp_frame(frame_size()), delayed_frame(run_length());

    vector<complexf> rotated_frame(run_length());

    noisy_sequence.resize(run_length());

    message_generation(message);
#ifdef DEBUG_GENERATOR_MD
    std::ofstream seq_log_file("message.txt", std::ios::app);
    for (const auto & val : message) {
        seq_log_file << val << "\n";
    }
    seq_log_file << endl;
    seq_log_file.close();
#endif

    nbldpc_encode(message, codeword);
#ifdef DEBUG_GENERATOR_MD
    seq_log_file.open("codeword.txt", std::ios::app);
    for (const auto & val : codeword) {
        seq_log_file << val << "\n";
    }
    seq_log_file << endl;
    seq_log_file.close();
#endif

    ccsk_modulate(codeword, ccsk_frame);
#ifdef DEBUG_GENERATOR_MD
    seq_log_file.open("ccsk_frame.txt", std::ios::app);
    for (const auto & val : ccsk_frame) {
        seq_log_file << val << "\n";
    }
    seq_log_file << endl;
    seq_log_file.close();
#endif

    overmodulate(ccsk_frame, qcsp_frame);
#ifdef DEBUG_GENERATOR_MD
    seq_log_file.open("qcsp_frame.txt", std::ios::app);
    for (const auto & val : qcsp_frame) {
        seq_log_file << val << "\n";
    }
    seq_log_file << endl;
    seq_log_file.close();
#endif

    add_delay(qcsp_frame, run_parameters, delayed_frame);
#ifdef DEBUG_GENERATOR_MD
    seq_log_file.open("delayed_frame.txt", std::ios::app | std::ios::out);
    for (const auto & val : delayed_frame) {
        seq_log_file << val << "\n";
    }
    seq_log_file << endl;
    seq_log_file.close();
#endif

    add_bpsk_phase_rotation(delayed_frame, run_parameters, rotated_frame);
#ifdef DEBUG_GENERATOR_MD
    seq_log_file.open("rotated_frame.txt", std::ios::app);
    for (const auto & val : rotated_frame) {
        seq_log_file << val.real() << " " << val.imag() << "\n";
    }
    seq_log_file << endl;
    seq_log_file.close();
#endif

    add_noise(rotated_frame, noisy_sequence);
#ifdef DEBUG_GENERATOR_MD
    seq_log_file.open("noisy_sequence.txt", std::ios::app | std::ios::out);
    for (const auto & val : noisy_sequence) {
        seq_log_file << val.real() << " " << val.imag() << "\n";
    }
    seq_log_file << endl;
    seq_log_file.close();
#endif

    return run_parameters;
}

template <>
CGenerator::frame_param_t CGenerator::generate<true>(vector<complexf> & noisy_sequence) {
    const frame_param_t run_parameters {0, 0, 0};

    vector<complexf> rotated_frame(run_length(), 0.f);
#if DEBUG_GENERATOR_FA
    std::ofstream seq_log_file("rotated_frame.txt", std::ios::app | std::ios::out);
    for (const auto & val : rotated_frame) {
        seq_log_file << val.real() << " " << val.imag() << "\n";
    }
    seq_log_file << endl;
    seq_log_file.close();
#endif

    noisy_sequence.resize(run_length());

    add_noise(rotated_frame, noisy_sequence);
#if DEBUG_GENERATOR_FA
    seq_log_file.open("noisy_sequence.txt", std::ios::app | std::ios::out);
    for (const auto & val : noisy_sequence) {
        seq_log_file << val.real() << " " << val.imag() << "\n";
    }
    seq_log_file << endl;
    seq_log_file.close();
#endif

    return run_parameters;
}
