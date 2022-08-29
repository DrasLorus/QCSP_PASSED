#include "./CGenerator.hpp"

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

        auto iter_dest = ccsk_frame.begin() + symbol * q();

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

        const auto iter_dest = qcsp_frame.begin() + (i << p());

        std::transform(iter_begin, iter_end, iter_dest, [over_symbol](bool b) { return !(b ^ over_symbol); });
    }
}

void CGenerator::add_delay(const vector<bool> & qcsp_frame, const frame_param_t & param, vector<bool> & delayed_frame) const {
    const int32_t delay = param.delay;

    std::fill(delayed_frame.begin(), delayed_frame.end(), 0);

    const auto iter_dest = delayed_frame.begin() + frame_size() * 2 + delay;

    std::copy(qcsp_frame.cbegin(), qcsp_frame.cend(), iter_dest);
} // Add delay to the frame

void CGenerator::add_bpsk_phase_rotation(const vector<bool> & delayed_frame, const frame_param_t & param, vector<complexf> rotated_frame) const {
    const float phase    = param.phase;
    const float rotation = param.rotation;

    for (int i = frame_size() * 2 + param.delay; i < frame_size() * 3 + param.delay; i++) {
        const float local_rotation = rotation * float(i) / float(q()) + phase;
        rotated_frame[i]           = std::exp(complexf(0, local_rotation)) * float(int32_t(delayed_frame[i]) * 2 - 1);
    }
} // Add BPSK + add phase and rotation if any

void CGenerator::add_noise(const vector<complexf> rotated_frame, vector<complexf> noisy_sequence) {
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
    nbldpc_encode(message, codeword);
    ccsk_modulate(codeword, ccsk_frame);
    overmodulate(ccsk_frame, qcsp_frame);

    add_delay(qcsp_frame, run_parameters, delayed_frame);
    add_bpsk_phase_rotation(delayed_frame, run_parameters, rotated_frame);
    add_noise(rotated_frame, noisy_sequence);

    return run_parameters;
}

template <>
CGenerator::frame_param_t CGenerator::generate<true>(vector<complexf> & noisy_sequence) {
    const frame_param_t run_parameters {0, 0, 0};

    vector<complexf> rotated_frame(run_length(), 0.f);

    noisy_sequence.resize(run_length());

    add_noise(rotated_frame, noisy_sequence);

    return run_parameters;
}
