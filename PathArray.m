% =========================================================================
%  PathArray.m  —  Free-path LUT for photon transport (Beer-Lambert)
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Look-up table of photon free-path lengths between successive scattering
%  events.  In a homogeneous medium with extinction coefficient mu_t,
%  Beer-Lambert gives the survival probability
%      P(no event up to s)  =  exp(-mu_t * s),
%  so the free-path PDF is the exponential
%      p(s)  =  mu_t * exp(-mu_t * s)  =  (1/mfpath) * exp(-s/mfpath)
%  where mfpath = 1/mu_t is the mean free path.  E[s] = mfpath.
%
%  The LUT is built over s in (0, 6*mfpath], which retains
%      1 - exp(-6) = 99.75 %  of the distribution.
%  Sampling is by uniform-integer index over the table (rejection-free).
%
%  ALTERNATIVE (no LUT): direct inverse-CDF sample  s = -mfpath * log(rand)
%  is mathematically exact in one line; the LUT trades memory for slightly
%  faster sampling and identical statistics.
%
%  References:
%    [1] Wang, L., Jacques, S.L., Zheng, L. (1995), "MCML—Monte Carlo
%        modeling of light transport in multi-layered tissues", Computer
%        Methods and Programs in Biomedicine 47, 131-146.  (Eq. for free
%        path sampling.)
%    [2] Ishimaru, A. (1978), "Wave Propagation and Scattering in Random
%        Media", Academic Press, Chap. 7.  (Beer-Lambert derivation.)
%    [3] Prahl, S.A. et al. (1989), SPIE Institute Series IS 5, 102-111.
% =========================================================================

function result = PathArray(step, mfpath)

    % Domain: 0 to 6 mfp (~99.75% of the exponential tail captured)
    s_max      = ceil(6 * mfpath);
    num        = ceil(s_max / step);
    stepp      = s_max / num;
    path_array = (stepp:stepp:s_max)';

    % Discrete exponential PDF  p(s) ~ exp(-s/mfpath)
    prob = exp(-path_array / mfpath);

    % Convert to integer bin occupancy; normalise so the smallest non-zero
    % bin keeps ~1 count (preserves the tail at finite LUT resolution)
    prob = ceil(1000 * prob / max(prob));

    % Cumulative -> LUT expansion
    Oprob  = cumsum(prob);
    totalN = Oprob(end);

    path_prob_array = zeros(1, totalN);
    path_prob_array(1:Oprob(1)) = path_array(1);
    for u = 2:num
        path_prob_array(Oprob(u-1)+1:Oprob(u)) = path_array(u);
    end

    result = path_prob_array;
end
