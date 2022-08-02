function filename = save_test_vectors(PN, best_N, rotation_span, p_omega, step_denominator, seed)
%SAVE_TEST_VECTORS Generate a set of test vectors and saves it in a matfile.
%   You must provide a PN sequence, a overmodulation sequence `best_N`, the rotation span (symbol 
%   rotation is in [-rotation_span/2, rotation_span/2[), the number of frequency hypothesis, the
%   step_denominator (maximum of error is pi/step_denominator) and a seed for the randomness. Be
%   aware that files grow with the length of PN and best_N, and with p_omega.
%
%   Usage:
%       fn = SAVE_TEST_VECTORS(PN, best_N, rotation_span, p_omega, step_denominator, seed)
%
%   Example:
%       fn = SAVE_TEST_VECTORS(PN, best_N, 2 * pi, 4, 4, 0) produce a file
%       "test_data_w4_step4_span0.5.mat" that contains:
%              cabs_max_sqr_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%              cabs_max_sqr_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%             cabs_max_sqr_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%             cabs_max_sqr_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%             cabs_max_sqrt_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%             cabs_max_sqrt_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%            cabs_max_sqrt_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%            cabs_max_sqrt_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%                   data_input_infdB_w4_q64_N60_1pi_2_n10: [192000x1 single]                                                                                               
%                   data_input_m10dB_w4_q64_N60_1pi_2_n30: [576000x1 single]                                                                                               
%                       deltas_infdB_w4_q64_N60_1pi_2_n10: [10x1     single]                                                                                               
%                       deltas_m10dB_w4_q64_N60_1pi_2_n30: [30x1     single]                                                                                               
%             iter_fcts_sqr_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%             iter_fcts_sqr_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%            iter_fcts_sqr_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%            iter_fcts_sqr_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%            iter_fcts_sqrt_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%            iter_fcts_sqrt_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%           iter_fcts_sqrt_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%           iter_fcts_sqrt_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%                 norms_sqr_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x1 single]                                                                                               
%                 norms_sqr_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x1 single]                                                                                               
%                norms_sqrt_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x1 single]                                                                                               
%                norms_sqrt_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x1 single]                                                                                               
%                    rotations_infdB_w4_q64_N60_1pi_2_n10: [10x1     single]                                                                                               
%                    rotations_m10dB_w4_q64_N60_1pi_2_n30: [30x1     single]                                                                                               
%                 score_sqr_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%                 score_sqr_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%                score_sqr_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%                score_sqr_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%                score_sqrt_l2_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%                score_sqrt_l2_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]                                                                                               
%               score_sqrt_raw_infdB_w4_q64_N60_1pi_2_n10: [192000x4 single]                                                                                               
%               score_sqrt_raw_m10dB_w4_q64_N60_1pi_2_n30: [576000x4 single]
%       and that can be loaded either using load or matfile (if supported).
%
%   See also LOAD, MATFILE

q = length(PN);
N = length(best_N);

rng(seed)

snr_tab  = [inf, -10];
snr_str  = {'infdB', 'm10dB'};
runs_tab = [10,   30];

    function stro = add_sqr(stri, squared)
        stro = cell(1, length(stri));
        for idx = 1 : length(stri)
            if squared
                stro{idx} = [stri{idx}, '_sqr'];
            else
                stro{idx} = [stri{idx}, '_sqrt'];
            end
        end
    end

    function stro = add_norm(stri, normed)
        stro = cell(1, length(stri));
        for idx = 1 : length(stri)
            if normed
                stro{idx} = [stri{idx}, '_l2'];
            else
                stro{idx} = [stri{idx}, '_raw'];
            end
        end
    end

    function stro = add_tail(stri, N, q, p_omega, rotation_span, i)
        stro = cell(1, length(stri));
        for idx = 1 : length(stri)
            [ratD, ratN] = rat(rotation_span / pi);
            stro{idx} = [stri{idx}, sprintf('_%s_w%i_q%i_N%i_%dpi_%d_n%i', ...
                snr_str{i}, p_omega, q, N, ratD, ratN, runs_tab(i))];
        end
    end



fieldroots = {'score', 'cabs_max', 'iter_fcts'};


fieldsin = {'data_input', 'deltas', 'rotations'};

%%%%%
%%
i = 1;

rng(seed)

[noisy_seq, ~, deltas, rotations, ~] = ...
    generate_noisy_sequence(runs_tab(i), PN, best_N, snr_tab(i), rotation_span);

local_fields = add_tail(fieldsin, N, q, p_omega, rotation_span, i);
local_data   = {single(noisy_seq), single(deltas), single(rotations)};

init_s = [local_fields; local_data];

struct_in_inf = struct(init_s{:});

%%
normed  = false;
squared = false;
[score, iter_fcts, cabs_max, ~] = ...
    compute_ts_score(p_omega, step_denominator, N, PN, normed, squared, noisy_seq);


local_fields = add_tail(add_norm(add_sqr(fieldroots, squared), normed), N, q, p_omega, rotation_span, i);
local_data   = {score, cabs_max, iter_fcts};

init_s = [local_fields; local_data];

struct_raw_sqrt_inf = struct(init_s{:});

%%
normed  = false;
squared = true;
[score, iter_fcts, cabs_max, ~] = ...
    compute_ts_score(p_omega, step_denominator, N, PN, normed, squared, noisy_seq);


local_fields = add_tail(add_norm(add_sqr(fieldroots, squared), normed), N, q, p_omega, rotation_span, i);
local_data   = {score, cabs_max, iter_fcts};

init_s = [local_fields; local_data];

struct_raw_sqr_inf = struct(init_s{:});

%%
normed  = true;
squared = false;
[score, iter_fcts, cabs_max, norms] = ...
    compute_ts_score(p_omega, step_denominator, N, PN, normed, squared, noisy_seq);


local_fields = add_tail(add_norm(add_sqr([fieldroots, {'norms'}], squared), normed), N, q, p_omega, rotation_span, i);
local_data   = {score, cabs_max, iter_fcts, norms};

init_s = [local_fields; local_data];

struct_l2_sqrt_inf = struct(init_s{:});

%%
normed  = true;
squared = true;
[score, iter_fcts, cabs_max, norms] = ...
    compute_ts_score(p_omega, step_denominator, N, PN, normed, squared, noisy_seq);


local_fields = add_tail(add_norm(add_sqr([fieldroots, {'norms'}], squared), normed), N, q, p_omega, rotation_span, i);
local_data   = {score, cabs_max, iter_fcts, norms};

init_s = [local_fields; local_data];

struct_l2_sqr_inf = struct(init_s{:});

%%
i = 2;

rng(seed)

[noisy_seq, ~, deltas, rotations, ~] = ...
    generate_noisy_sequence(runs_tab(i), PN, best_N, snr_tab(i), rotation_span);

local_fields = add_tail(fieldsin, N, q, p_omega, rotation_span, i);
local_data   = {noisy_seq, deltas, rotations};

init_s = [local_fields; local_data];

struct_in_m10 = struct(init_s{:});

%%
normed  = false;
squared = false;
[score, iter_fcts, cabs_max, ~] = ...
    compute_ts_score(p_omega, step_denominator, N, PN, normed, squared, noisy_seq);


local_fields = add_tail(add_norm(add_sqr(fieldroots, squared), normed), N, q, p_omega, rotation_span, i);
local_data   = {score, cabs_max, iter_fcts};

init_s = [local_fields; local_data];

struct_raw_sqrt_m10 = struct(init_s{:});

%%
normed  = false;
squared = true;
[score, iter_fcts, cabs_max, ~] = ...
    compute_ts_score(p_omega, step_denominator, N, PN, normed, squared, noisy_seq);


local_fields = add_tail(add_norm(add_sqr(fieldroots, squared), normed), N, q, p_omega, rotation_span, i);
local_data   = {score, cabs_max, iter_fcts};

init_s = [local_fields; local_data];

struct_raw_sqr_m10 = struct(init_s{:});

%%
normed  = true;
squared = false;
[score, iter_fcts, cabs_max, norms] = ...
    compute_ts_score(p_omega, step_denominator, N, PN, normed, squared, noisy_seq);


local_fields = add_tail(add_norm(add_sqr([fieldroots, {'norms'}], squared), normed), N, q, p_omega, rotation_span, i);
local_data   = {score, cabs_max, iter_fcts, norms};

init_s = [local_fields; local_data];

struct_l2_sqrt_m10 = struct(init_s{:});

%%
normed  = true;
squared = true;
[score, iter_fcts, cabs_max, norms] = ...
    compute_ts_score(p_omega, step_denominator, N, PN, normed, squared, noisy_seq);


local_fields = add_tail(add_norm(add_sqr([fieldroots, {'norms'}], squared), normed), N, q, p_omega, rotation_span, i);
local_data   = {score, cabs_max, iter_fcts, norms};

init_s = [local_fields; local_data];

struct_l2_sqr_m10 = struct(init_s{:});

full_struct = combineStructs(...
    combineStructs(...
        combineStructs(...
            combineStructs(...
                combineStructs(...
                    combineStructs(...
                        combineStructs(...
                            combineStructs(...
                                combineStructs(struct_in_inf, struct_in_m10), ...
                                struct_raw_sqrt_inf), ...
                            struct_raw_sqrt_m10), ...
                        struct_raw_sqr_inf), ...
                    struct_raw_sqr_m10), ...
                struct_l2_sqrt_inf), ...
            struct_l2_sqrt_m10), ...
        struct_l2_sqr_inf), ...
    struct_l2_sqr_m10);

filename = sprintf("test_data_q%i_N%i_w%i_step%i_span%3.1f.mat", q, N, p_omega, step_denominator, rotation_span / pi);

if exist('OCTAVE_VERSION', 'builtin')
    save(filename, '-v7', '-struct', 'full_struct')
else
    save(filename, '-v7.3', '-struct', 'full_struct')
end

end
