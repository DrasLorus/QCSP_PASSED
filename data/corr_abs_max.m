function cabs_max = corr_abs_max(iterative_factors, PN, norm_factors, normed, squared)
%  CORR_ABS_MAX Correlate with PN following the time sliding method and output the absolute maximum
%   It consumes an iterative factor column vector (iterative_factors) as produced be the time 
%   sliding method and correlate it independently for each frequency hypothesis with the PN
%   sequence. The correlation is iterative, thus to use a different PN or work on uncorrelated data,
%   "clear corr_abs_max" must be run to clear the persistent correlation register vector, the
%   behavior is undefined otherwise. The absolute maxima of correlation (one for each frequency) can
%   then be optionaly normalize br norm_factors and/or squared.
%
%   Usage:
%       cabs_max = corr_abs_max(iterative_factors, PN, norm_factors, normed, squared)
%         outputs the absolute maximum of correlation, normalized by norm_factors if normed is True,
%         and squared if squared is True.
%
%   Critical Notice:
%       This function may be a bit technical to use out of context, due to its iterative nature. It
%       is better to use it from a dedicated function like compute_ts_score. To compute a one shot
%       correlation between a vector X and PN, it is advised to use the FFT method, e.g.
%           cabs_max = max(abs(ifft(fft(X, [], 2) .* conj(fft(PN64, [], 2)), [], 2)), [], 2)
%
%   See also GENERATE_NOISY_SEQUENCE, COMPUTE_TS_SCORE, FFT, IFFT, MAX

persistent corr_i_r
if isempty(corr_i_r)
    corr_i_r = zeros(length(iterative_factors), length(PN));
end

corr_ip1 = corr_i_r + PN .* iterative_factors;

corr_i_r = circshift(corr_ip1, 1 ,2);


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
