% =========================================================================
%  Mie_abcd.m  —  Mie expansion coefficients a_n, b_n, c_n, d_n
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Returns a 4 x nmax matrix [a_n; b_n; c_n; d_n] of complex Mie
%  coefficients for relative refractive index m = m' + i*m" and size
%  parameter x = 2*pi*a*n_medium/lambda.
%
%  Equations  (Bohren & Huffman 1983, p.100):
%      a_n = ( m^2 psi_n(mx) psi_n'(x) - psi_n(x) psi_n'(mx) ) /
%            ( m^2 psi_n(mx) xi_n'(x) -  xi_n(x) psi_n'(mx) )
%      b_n = (     psi_n(mx) psi_n'(x) - psi_n(x) psi_n'(mx) ) /
%            (     psi_n(mx) xi_n'(x)  -  xi_n(x) psi_n'(mx) )
%      c_n,d_n: internal-field coefficients (B&H 4.50, 4.51)
%
%  where psi_n, xi_n are Riccati-Bessel functions built from spherical
%  Bessel/Hankel functions via
%      psi_n(z) = z j_n(z),   xi_n(z) = z h_n^(1)(z) = psi_n(z) + i chi_n(z).
%
%  Truncation:  nmax = round(2 + x + 4*x^(1/3))  (Wiscombe 1980).
%
%  References:
%    [1] Bohren, C.F. & Huffman, D.R. (1983), "Absorption and Scattering of
%        Light by Small Particles", Wiley, pp. 100, 477.
%    [2] Wiscombe, W.J. (1980), "Improved Mie scattering algorithms",
%        Appl. Opt. 19, 1505-1509  (nmax truncation criterion).
%    [3] Maetzler, C. (2002), "MATLAB Functions for Mie Scattering and
%        Absorption", IAP Research Report 2002-08, University of Bern.
%
%  Original code: Christian Maetzler, June 2002.  Modified: Dongna, Apr 2012.
% =========================================================================

function result = Mie_abcd(m, x)

    nmax = round(2 + x + 4*x^(1/3));       % Wiscombe 1980
    n    = (1:nmax);  nu = n + 0.5;
    z    = m .* x;   m2 = m .* m;

    % Riccati-Bessel functions via half-integer Bessel functions:
    %   j_n(z) = sqrt(pi/(2z)) * J_{n+1/2}(z)
    sqx = sqrt(0.5*pi ./ x);
    sqz = sqrt(0.5*pi ./ z);

    bx  = besselj(nu, x) .* sqx;            % psi_n(x)/x  (order n)
    bz  = besselj(nu, z) .* sqz;            % psi_n(mx)/mx
    yx  = bessely(nu, x) .* sqx;            % chi_n(x)/x
    hx  = bx + 1i*yx;                        % xi_n(x)/x  =  h_n^(1)(x)

    % displaced (n-1) arrays for the derivative recursion psi_n' = ...
    b1x = [sin(x)/x,  bx(1:nmax-1)];
    b1z = [sin(z)/z,  bz(1:nmax-1)];
    y1x = [-cos(x)/x, yx(1:nmax-1)];
    h1x = b1x + 1i*y1x;

    % Riccati-Bessel derivatives (B&H eq. 4.89)
    ax  = x .* b1x - n .* bx;
    az  = z .* b1z - n .* bz;
    ahx = x .* h1x - n .* hx;

    % Mie coefficients  (B&H eqs. 4.56, 4.57, 4.50, 4.51)
    an = (m2 .* bz .* ax  - bx  .* az) ./ (m2 .* bz .* ahx - hx .* az);
    bn = (      bz .* ax  - bx  .* az) ./ (      bz .* ahx - hx .* az);
    cn = (      bx .* ahx - hx  .* ax) ./ (      bz .* ahx - hx .* az);
    dn = m .* (bx .* ahx - hx  .* ax) ./ (m2 .* bz .* ahx - hx .* az);

    result = [an; bn; cn; dn];
end
