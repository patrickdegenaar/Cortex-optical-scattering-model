% =========================================================================
%  Mie.m  —  Mie scattering / extinction / absorption efficiencies
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Computes, for relative refractive index m = m' + i*m" and size parameter
%  x = 2*pi*a*n_medium/lambda,
%      Q_ext, Q_sca, Q_abs, Q_back, g (= <cos theta>), Q_back/Q_sca
%  using the complex Mie coefficients a_n, b_n returned by Mie_abcd.
%
%  Result vector (1x9):
%      [m'   m"   x   Q_ext   Q_sca   Q_abs   Q_back   g   Q_back/Q_sca]
%
%  Equations (Bohren & Huffman 1983, "Absorption and Scattering of Light by
%  Small Particles", Wiley):
%    Q_ext   = (2/x^2) * sum_n (2n+1) * Re(a_n + b_n)            (B&H 4.61)
%    Q_sca   = (2/x^2) * sum_n (2n+1) * (|a_n|^2 + |b_n|^2)      (B&H 4.62)
%    Q_abs   = Q_ext - Q_sca
%    Q_back  = (1/x^2) * |sum_n (2n+1) * (-1)^n * (a_n - b_n)|^2 (B&H 4.83)
%    g       = (4/(x^2 Q_sca)) *
%              { sum_n  n(n+2)/(n+1) * Re(a_n a*_{n+1} + b_n b*_{n+1})
%               + sum_n (2n+1)/(n(n+1)) * Re(a_n b*_n) }         (B&H 4.64)
%
%  References:
%    [1] Bohren, C.F. & Huffman, D.R. (1983), "Absorption and Scattering of
%        Light by Small Particles", Wiley, pp. 103, 119-122, 477.
%    [2] Maetzler, C. (2002), "MATLAB Functions for Mie Scattering and
%        Absorption", IAP Research Report 2002-08, University of Bern.
%    [3] van de Hulst, H.C. (1957), "Light Scattering by Small Particles",
%        Wiley (Dover reprint 1981).
%
%  Original code: Christian Maetzler, May 2002.  Modified: Dongna, Apr 2012.
% =========================================================================

function result = Mie(m, x)

    if x == 0                       % singularity at x = 0
        result = [real(m) imag(m) 0 0 0 0 0 0 1.5];
    elseif x > 0                    % normal case
        nmax = round(2 + x + 4*x^(1/3));        % Wiscombe (1980) criterion
        n1   = nmax - 1;
        n    = 1:nmax;
        cn   = 2*n + 1;
        c1n  = n.*(n+2)./(n+1);                  % weights for asym sum
        c2n  = cn ./ n ./ (n+1);
        x2   = x*x;

        f    = Mie_abcd(m, x);                   % Mie coeffs a_n,b_n,c_n,d_n
        anp  = real(f(1,:));  anpp = imag(f(1,:));
        bnp  = real(f(2,:));  bnpp = imag(f(2,:));

        % displaced (n -> n+1) coefficients used in the g sum (B&H p.120)
        g1(1:4,nmax) = [0; 0; 0; 0];
        g1(1,1:n1)   = anp(2:nmax);
        g1(2,1:n1)   = anpp(2:nmax);
        g1(3,1:n1)   = bnp(2:nmax);
        g1(4,1:n1)   = bnpp(2:nmax);

        % --- Q_ext  (B&H 4.61) ---
        dn   = cn .* (anp + bnp);
        qext = 2*sum(dn) / x2;

        % --- Q_sca  (B&H 4.62) ---
        en   = cn .* (anp.*anp + anpp.*anpp + bnp.*bnp + bnpp.*bnpp);
        qsca = 2*sum(en) / x2;

        qabs = qext - qsca;

        % --- Q_back (B&H 4.83) ---
        fn   = (f(1,:) - f(2,:)) .* cn;
        gn   = (-1).^n;
        f(3,:) = fn .* gn;
        qb   = (sum(f(3,:)) * sum(f(3,:))') / x2;

        % --- asymmetry parameter g = <cos theta>  (B&H 4.64) ---
        asy1 = c1n .* (anp.*g1(1,:) + anpp.*g1(2,:) + bnp.*g1(3,:) + bnpp.*g1(4,:));
        asy2 = c2n .* (anp.*bnp + anpp.*bnpp);
        asy  = 4/x2 * sum(asy1 + asy2) / qsca;

        qratio = qb / qsca;
        result = [real(m) imag(m) x qext qsca qabs qb asy qratio];
    end
end
