function save_test_vectors(PN, best_N, rotation_span, p_omega, step_denominator, seed)
%SAVE_TEST_VECTORS Generate aset of test vectors and saves it in an HDF5 matfile.
%   You must provide a PN sequence, a overmodulation sequence `best_N`, the rotation span (symbol 
%   rotation is in [-rotation_span/2, rotation_span/2[), the number of frequency hypothesis, the
%   step_denominator (maximum of error is pi/step_denominator) and a seed for the randomness. Be
%   aware that files grow with the length of PN and best_N, and with p_omega.

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
            stro{idx} = [stri{idx}, sprintf('_%s_w%i_q%i_N%d_%dpi_n%i', ...
                snr_str{i}, p_omega, q, N, rotation_span / pi, runs_tab(i))];
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
local_data   = {noisy_seq, deltas, rotations};

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

save(sprintf("test_data_w%i_pi_%i.mat", p_omega, step_denominator), '-v7.3', '-struct', 'full_struct')

end

