% =========================================================================
%  Parameters.m  —  Derived Mie inputs and per-particle cross-section
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Returns the canonical Mie inputs for a sphere of diameter d (um) in a
%  medium of refractive index n_e at wavelength lambda (um):
%
%      m   = n_p / n_e                                  (relative index)
%      x   = 2 pi a n_e / lambda  =  pi d / (lambda / n_e)   (size parameter)
%      pcs = pi d^2 / 4           (geometric cross-section, um^2)
%      pV  = (4/3) pi (d/2)^3     (sphere volume, um^3)
%
%  m may be complex (m = m' + i*m");  the imaginary part is the absorption
%  index of the scatterer material.
%
%  Bulk coefficients are then formed in the caller as
%      mu_ext = rho * pcs * Q_ext(m,x),
%      mu_abs = rho * pcs * Q_abs(m,x),
%      mu_sca = mu_ext - mu_abs,
%  with rho the number density of scatterers.
%
%  References:
%    [1] Bohren, C.F. & Huffman, D.R. (1983), "Absorption and Scattering of
%        Light by Small Particles", Wiley, Sec. 4.1.
%    [2] Jacques, S.L. (2013), "Optical properties of biological tissues:
%        a review", Phys. Med. Biol. 58, R37-R61.
% =========================================================================

function [m, x, pcs, pV] = Parameters(lamda, d, ne, np)

    pcs = pi*(d^2)/4;
    pV  = (4*pi*(d/2)^3)/3;

    m = np / ne;
    x = (pi*d) / (lamda/ne);
end
