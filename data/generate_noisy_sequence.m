function [noisy_seq, frames, deltas, rotations, codewords] = generate_noisy_sequence(runs, pn, ovmod, snr, rotation_span)

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
