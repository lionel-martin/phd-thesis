function [ filt_sig ] = gsp_filter_new_ideal_lowpass( G, thres, signal, param )
%   GSP_FILTER_NEW_IDEAL_LOWPASS filters a signal with the approximation of
%   the ideal lowpass filter using the approximation of the sign function
%   proposed by Allen-Zhu et al.
%
%   filt_sig = gsp_filter_new_ideal_lowpass(G, thres, signal)
%   generates the response of a signal with a lowpass filter of threshold
%   thres.
%
%   filt_sig = gsp_filter_new_ideal_lowpass(G, thres, signal, param)
%   allows to pass additional parameters:
%   * 'order' handles the order of the polynomial approximation
%   * 'kap' is the kappa parameter of the Allen-Zhu approximation.
%
%   Developed under Matlab version 8.5.0.197613 (R2015a)
%   Created by Lionel Martin

    assert(size(signal, 1) == G.N);
    n = G.N;

    if nargin < 4, param = struct(); end
    if ~isfield(param, 'kap'), param.kap = 0.01; end
    if ~isfield(param, 'order'), param.order = 50; end

    if ~isfield(G, 'lmax'), G = gsp_estimate_lmax(G); end

    tau = max(thres, G.lmax-thres);
    shift = (thres * speye(n) - G.L);

    h = @(x) 1./sqrt((1 + param.kap - x) / 2); 
    c = gsp_cheby_coeff([-1, 1], h, param.order);

    T_old = shift * signal;
    T_cur = (1 + param.kap) * T_old - 2 / tau^2 * shift * (shift * T_old);

    filt_sig = signal / 2 + (c(1) / (4 * tau)) * T_old; % degree 0 term
    filt_sig = filt_sig + (c(2) / (2 * tau)) * T_cur; % degree 1 term

    for ii = 3:param.order+1
        T_new = (2 + 2 * param.kap) * T_cur - 4 / tau^2 * shift * (shift * T_cur) - T_old;
        filt_sig = filt_sig + c(ii)/(2*tau) * T_new;
        T_old = T_cur;
        T_cur = T_new;
    end
end
