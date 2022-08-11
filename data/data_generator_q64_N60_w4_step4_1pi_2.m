
load('parameters_20210903.mat', 'PN64', 'best_N');
fnm = save_test_vectors(PN64, best_N, pi / 2, 4, 4, 0);
fprintf('%s generated.\n', fnm)

