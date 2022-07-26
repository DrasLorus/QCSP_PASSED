function [x_idxs, y_idxs, lin_indeces] = get_max_indeces(score, runs, N, q)

run_length = size(score, 1) / runs;
p_omega = size(score, 2);

maxv = zeros(runs, p_omega);
idxs = zeros(runs, p_omega);

for j = 1 : runs
    [maxv(j, :), idxs(j, :)] = max(score((j - 1) * run_length + 1 : j * run_length, :));
end
[~, I] = max(maxv, [], 2);

x_idxs = idxs(sub2ind(size(idxs), (1 : runs)', I(:))) + (0 : runs - 1)' * run_length;
y_idxs = I;
lin_indeces = sub2ind(size(score), x_idxs, y_idxs);
