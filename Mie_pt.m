% =========================================================================
%  Mie_pt.m  —  Angular functions pi_n(cos theta) and tau_n(cos theta)
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Computes the angular functions used in the Mie scattering amplitudes
%  S_1, S_2.  For u = cos(theta) and integer n >= 1:
%      pi_n(u)  = P_n^1(u) / sin(theta)
%      tau_n(u) = d/d(theta) [P_n^1(u)]
%  where P_n^1 is the associated Legendre function of degree n, order 1.
%
%  Computed by recurrence  (Bohren & Huffman 1983, eqs. 4.47):
%      pi_1   = 1,         tau_1 = u
%      pi_2   = 3u,        tau_2 = 3 cos(2 theta) = 3(2u^2 - 1)
%      pi_n   = (2n-1)/(n-1) * u * pi_{n-1}   -   n/(n-1) * pi_{n-2}
%      tau_n = n*u*pi_n  -  (n+1)*pi_{n-1}
%
%  References:
%    [1] Bohren, C.F. & Huffman, D.R. (1983), "Absorption and Scattering of
%        Light by Small Particles", Wiley, pp. 94-95.
%    [2] Maetzler, C. (2002), "MATLAB Functions for Mie Scattering and
%        Absorption", IAP Research Report 2002-08, University of Bern.
% =========================================================================

function result = Mie_pt(u, nmax)

    p    = zeros(1, nmax);
    t    = zeros(1, nmax);

    p(1) = 1;
    t(1) = u;
    p(2) = 3*u;
    t(2) = 3*cos(2*acos(u));                  % == 3*(2u^2 - 1)

    for n1 = 3:nmax
        p(n1) = (2*n1 - 1)/(n1 - 1) * p(n1-1) * u  -  n1/(n1 - 1) * p(n1-2);
        t(n1) = n1 * u * p(n1)                     -  (n1 + 1) * p(n1-1);
    end

    result = [p; t];
end
