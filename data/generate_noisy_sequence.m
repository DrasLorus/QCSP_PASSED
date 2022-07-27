function [noisy_seq, frames, deltas, rotations, codewords] = generate_noisy_sequence(runs, pn, ovmod, snr, rotation_span)
%  GENERATE_NOISY_SEQUENCE Generate a random noisy CCSK sequence
%   The number of runs (runs) must be provided, with a PN sequence (pn), an overmodulation sequence
%   (ovmod), the targeted SNR (snr) and the rotation span (rotation_span).
%   The sequence represent runs CCSK frames of length(ovmod) symbols of length(pn) chips separated
%   by at most 4 * length(ovmod) * length(pn) zeros and at least 3 * length(ovmod) * length(pn)
%   zeros. The begining index (the time inacurracy, or Delta) is random. Each sequence is randomly
%   rotated bsuch that a symbol can rotate from -rotation_span / 2 to rotation_span / 2. Finally, a 
%   complex white gaussian noise is added, with a power ensuring an Signal-to-Noise Ration of snr.
%
%   Usage:
%       noisy_seq =  generate_noisy_sequence(runs, pn, ovmod, snr, rotation_span)
%         outputs only the noisy sequence
%
%       [noisy_seq, frames, deltas, rotations, codewords] = ...
%           generate_noisy_sequence(runs, pn, ovmod, snr, rotation_span)
%         also outputs the clear frames, deltas, rotations and sent codewords.
%
%   See also CORR_ABS_MAX, COMPUTE_TS_SCORE


q = length(pn);
N = length(ovmod);
% p_theta = 5;

% snr = 10;

if size(pn, 1) == 1
    pn = reshape(pn, [], 1);
end

sigma   = sqrt(10^(-snr/10));
sigma_c = sigma / sqrt(2);

frame_length = N * q;

codewords = randi([0, q - 1], runs, N);
frames    = zeros(frame_length, runs);
for run = 1 : runs
    for symbol = 1 : N
        frames((symbol - 1) * q + 1 : symbol * q, run) = circshift(pn, -codewords(run, symbol)) .* ovmod(symbol);
    end
end

symbol_rot = -rotation_span + 2 * rotation_span * rand(1, runs);

delta = randi([0, q - 1], 1, runs);

deltas    = single(reshape(delta, runs, 1));
rotations = single(reshape(symbol_rot, runs, 1));

run_length = 5 * N * q;
seq        = zeros(run_length * runs, 1);
for n_run = 1 : runs
    idx_seq = (n_run - 1) * run_length + 1 : n_run * run_length;
    seq(idx_seq) = [zeros(2 * N * q + delta(n_run), 1) ; frames(:, n_run) .* exp(symbol_rot(n_run) * 1i / q *...
        reshape(1 : frame_length, frame_length, 1)) ; zeros(2 * N * q - delta(n_run), 1)];
    % seq = seq
    %seq = seq + sigma * (randn(1, length(seq)));
end

noisy_seq = single(seq + sigma_c * (randn(size(seq)) + 1i * randn(size(seq))));
