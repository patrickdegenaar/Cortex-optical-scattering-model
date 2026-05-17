% =========================================================================
%  Mie_scat_angle_array.m  —  Scattering-angle LUT from the Mie phase fn
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Builds a look-up table of polar scattering angles theta in [0, pi]
%  sampled in proportion to the Mie polar-angle PDF
%      w(theta)  =  p(cos theta) * sin(theta),
%  i.e. the phase function returned by Mie_scat_pattern multiplied by the
%  solid-angle (sin theta) Jacobian.  The azimuthal angle phi is uniformly
%  sampled in [0, 2*pi).
%
%  WHY THE sin(theta) JACOBIAN?
%  The solid-angle element is  d Omega = sin(theta) d theta d phi, so a
%  density per steradian, p(cos theta), maps to a per-theta density via
%      p_theta(theta) = p(cos theta) * sin(theta)
%  after integrating over phi.  Omitting sin(theta) over-weights the polar
%  caps and gives wrong scattering statistics.
%
%  References for sin(theta) Jacobian + Mie sampling in MC photon transport:
%    [1] Wang, L., Jacques, S.L., Zheng, L. (1995), "MCML—Monte Carlo
%        modeling of light transport in multi-layered tissues", Computer
%        Methods and Programs in Biomedicine 47, 131-146.
%    [2] Prahl, S.A., Keijzer, M., Jacques, S.L., Welch, A.J. (1989),
%        "A Monte Carlo Model of Light Propagation in Tissue", SPIE
%        Institute Series IS 5, 102-111.
%    [3] Bohren, C.F. & Huffman, D.R. (1983), Sec. 4.7.
% =========================================================================

function result = Mie_scat_angle_array(m, a)

    TH    = 0:pi/180:pi;                  % 1-degree resolution in theta
    nTH   = length(TH);
    Mie_p = zeros(1, nTH);

    for u = 1:nTH
        Mie_p(u) = Mie_scat_pattern(m, a, cos(TH(u)));
    end

    % --- sin(theta) solid-angle Jacobian ---
    Mie_p = Mie_p .* abs(sin(TH));

    prob = Mie_p ./ sum(Mie_p);
    prob = ceil(60000 .* prob);            % integer bin occupancy

    Oprob  = cumsum(prob);
    totalN = Oprob(end);

    scat_angle_array = zeros(1, totalN);
    scat_angle_array(1:Oprob(1)) = TH(1);
    for u = 2:nTH
        scat_angle_array(Oprob(u-1)+1:Oprob(u)) = TH(u);
    end

    % Azimuth phi: uniform on [0, 2*pi)
    azimu_angle_array = (1:totalN) * 2*pi / totalN;

    result = [scat_angle_array; azimu_angle_array];
end
