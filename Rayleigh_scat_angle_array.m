% =========================================================================
%  Rayleigh_scat_angle_array.m  —  Rayleigh scattering-angle LUT
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Look-up table of scattering angles drawn from the unpolarised Rayleigh
%  phase function
%      p(cos theta)  proportional to  (1 + cos^2(theta))
%  with the per-theta sampling weight (including the sin(theta) solid-angle
%  Jacobian)
%      w(theta)  =  (1 + cos^2(theta)) * sin(theta).
%  The azimuthal angle phi is uniform on [0, 2*pi).
%
%  The Rayleigh regime applies for size parameter x = 2*pi*a/lambda << 1
%  (and m*x << 1).  The phase function above is the dipole pattern, with
%  forward / backward maxima of  1 + 1 = 2  and a side minimum of  1.
%
%  References:
%    [1] Lord Rayleigh / J.W. Strutt (1871), "On the light from the sky,
%        its polarization and colour", Phil. Mag. 41, 107-120, 274-279.
%    [2] Bohren, C.F. & Huffman, D.R. (1983), "Absorption and Scattering of
%        Light by Small Particles", Wiley, p. 132.
%    [3] van de Hulst, H.C. (1957), Sec. 6.3.
% =========================================================================

function result = Rayleigh_scat_angle_array(Angle_step)

    Arc_step  = Angle_step * pi/180;
    num_angle = ceil(pi/Arc_step);
    TH        = 0:(pi/num_angle):pi;

    % Rayleigh weight  =  (1 + cos^2 theta) * sin theta
    prob = (1 + cos(TH).^2) .* abs(sin(TH));
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
