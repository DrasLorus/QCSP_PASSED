function [noisy_seq, frames, deltas, rotations, codewords] = generate_noisy_sequence(runs, pn, ovmod, snr, rotation_span)

q = length(pn);
N = length(ovmod);
% p_theta = 5;

% snr = 10;

sigma   = sqrt(10^(-snr/10));
sigma_c = sigma / sqrt(2);


codewords = randi([0, q - 1], runs, N);
frames    = zeros(runs, N * q);
for run = 1 : runs
    for symbol = 1 : N
        frames(run, (symbol - 1) * q + 1 : symbol * q) = circshift(pn, -codewords(run, symbol)) .* ovmod(symbol);
    end
end

symbol_rot = -rotation_span + 2 * rotation_span * rand(1, runs);

delta = randi([0, q - 1], 1, runs);

deltas    = delta;
rotations = symbol_rot;

run_length = 5 * N * q;
seq        = zeros(1,  run_length * runs);
for n_run = 1 : runs
    idx_seq = (n_run - 1) * run_length + 1 : n_run * run_length;
    seq(idx_seq) = [ zeros(1, 2 * N * q + delta(n_run)) frames(n_run , :).*exp(symbol_rot(n_run) * 1i / q *...
        (1:size(frames, 2))) zeros(1, 2 * N * q - delta(n_run))];
    % seq = seq
    %seq = seq + sigma * (randn(1, length(seq)));
end

noisy_seq = seq + sigma_c * (randn(1, length(seq)) + 1i * randn(1, length(seq)));
