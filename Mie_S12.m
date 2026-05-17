% =========================================================================
%  Mie_S12.m  —  Complex scattering amplitudes S_1(theta), S_2(theta)
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  For a sphere of size parameter x and relative refractive index m, and
%  for u = cos(theta), the far-field scattering amplitudes are
%      S_1(theta) = sum_n  (2n+1)/(n(n+1))  *  [ a_n pi_n(u) + b_n tau_n(u) ]
%      S_2(theta) = sum_n  (2n+1)/(n(n+1))  *  [ a_n tau_n(u) + b_n pi_n(u) ]
%  (Bohren & Huffman 1983, eqs. 4.74).
%
%  These are the perpendicular (S_1) and parallel (S_2) components of the
%  scattered electric field, normalised so that the differential scattering
%  cross-section per polarisation is
%      d sigma / d Omega = |S|^2 / k^2.
%
%  References:
%    [1] Bohren, C.F. & Huffman, D.R. (1983), pp. 111-114.
%    [2] Maetzler, C. (2002), IAP Research Report 2002-08, Univ. Bern.
%
%  Original code: Christian Maetzler, May 2002.  Modified: Dongna, Apr 2012.
% =========================================================================

function result = Mie_S12(m, x, u)

    nmax = round(2 + x + 4*x^(1/3));

    abcd = Mie_abcd(m, x);
    an   = abcd(1,:);
    bn   = abcd(2,:);

    pt   = Mie_pt(u, nmax);
    pin  = pt(1,:);
    tin  = pt(2,:);

    n    = 1:nmax;
    n2   = (2*n + 1) ./ (n .* (n + 1));      % series weight (B&H eq.4.74)
    pin  = n2 .* pin;
    tin  = n2 .* tin;

    S1   = an*pin' + bn*tin';                 % (B&H eq. 4.74a)
    S2   = an*tin' + bn*pin';                 % (B&H eq. 4.74b)

    result = [S1; S2];
end
