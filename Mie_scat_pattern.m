% =========================================================================
%  Mie_scat_pattern.m  —  Normalised Mie phase function p(cos theta)
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Returns the unit-integral phase function (per steradian)
%      p(cos theta) =  ( |S_1|^2 + |S_2|^2 )  /  ( 2  k^2  C_sca )
%  using
%      k^2 * C_sca  =  2 pi * sum_n (2n+1) (|a_n|^2 + |b_n|^2)
%  so that
%      p(cos theta) =  ( |S_1|^2 + |S_2|^2 )
%                     / ( 4 pi  *  sum_n (2n+1) (|a_n|^2 + |b_n|^2) ).
%
%  Normalisation:  integral over the full sphere == 1, i.e.
%      integral_0^pi p(cos theta) * 2 pi sin(theta) d theta = 1.
%
%  IMPORTANT (callers): p(cos theta) is a per-steradian density.  To sample
%  the polar angle theta, one must include the sin(theta) Jacobian:
%      w(theta) = p(cos theta) * sin(theta).
%  See Mie_scat_angle_array.m for that step.
%
%  References:
%    [1] Bohren, C.F. & Huffman, D.R. (1983), "Absorption and Scattering of
%        Light by Small Particles", Wiley, pp. 111-114, 381.
%    [2] van de Hulst, H.C. (1957), Sec. 9.32.
%    [3] Wang, L., Jacques, S.L. & Zheng, L. (1995), "MCML — Monte Carlo
%        modeling of light transport in multi-layered tissues", Computer
%        Methods and Programs in Biomedicine 47, 131-146.  (Phase-function
%        use in tissue Monte Carlo.)
% =========================================================================

function result = Mie_scat_pattern(m, x, u)

    nmax = round(2 + x + 4*x^(1/3));

    S12  = Mie_S12(m, x, u);
    S1   = S12(1,:);
    S2   = S12(2,:);

    abcd = Mie_abcd(m, x);
    an   = abcd(1,:);
    bn   = abcd(2,:);

    numerator   = abs(S1)^2 + abs(S2)^2;
    dnmntr      = (2*(1:nmax) + 1) .* ( abs(an(1:nmax)).^2 + abs(bn(1:nmax)).^2 );
    denominator = 4*pi * sum(dnmntr);

    result = numerator / denominator;
end
