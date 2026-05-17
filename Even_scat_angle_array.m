% =========================================================================
%  Even_scat_angle_array.m  —  Isotropic scattering-angle LUT
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Builds a look-up table for ISOTROPIC scattering (uniform on the sphere).
%  Despite the name "Even", the polar-angle distribution is NOT uniform in
%  theta — uniform-on-sphere implies
%      p(theta) = (1/2) * sin(theta),   theta in [0, pi]
%  because the solid-angle element  d Omega = sin(theta) d theta d phi.
%  The azimuthal angle phi is uniform on [0, 2*pi).
%
%  This is the Henyey-Greenstein phase function in the g -> 0 limit.
%
%  References:
%    [1] Wang, L., Jacques, S.L., Zheng, L. (1995), "MCML—Monte Carlo
%        modeling of light transport in multi-layered tissues", Computer
%        Methods and Programs in Biomedicine 47, 131-146.
%    [2] Ishimaru, A. (1978), "Wave Propagation and Scattering in Random
%        Media", Academic Press, Sec. 7-1 (isotropic scattering).
% =========================================================================

function result = Even_scat_angle_array(Angle_step)

    Arc_step  = Angle_step * pi/180;
    num_angle = ceil(pi/Arc_step);
    TH        = 0:(pi/num_angle):pi;

    % --- sin(theta) Jacobian for uniform-on-sphere sampling ---
    prob = abs(sin(TH));
    prob = prob ./ sum(prob);
    prob = ceil(40000 .* prob);

    Oprob  = cumsum(prob);
    totalN = Oprob(end);

    theta_array = zeros(1, totalN);
    theta_array(1:Oprob(1)) = TH(1);
    for u = 2:length(TH)
        theta_array(Oprob(u-1)+1:Oprob(u)) = TH(u);
    end

    % Azimuth phi: uniform on [0, 2*pi)
    phai_array = (1:totalN) * 2*pi / totalN;

    result = [theta_array; phai_array];
end
