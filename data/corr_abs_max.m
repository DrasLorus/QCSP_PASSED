function cabs_max = corr_abs_max(iterative_factors, PN, norm_factors, normed, squared)

persistent corr_i_r
if isempty(corr_i_r)
    corr_i_r = zeros(length(iterative_factors), length(PN));
end

corr_ip1 = corr_i_r + PN .* iterative_factors;

if exist('OCTAVE_VERSION', 'builtin')
	corr_i_r = correlation_shift_oct(corr_ip1);
else
	corr_i_r = circshift(corr_ip1, 1 ,2);
end

if squared
    vec_new_max = max(corr_ip1 .* conj(corr_ip1), [], 2);
else
    vec_new_max = max(abs(corr_ip1), [], 2);
end

if normed
   vec_new_max = vec_new_max .* norm_factors;
end

cabs_max = vec_new_max;

end
