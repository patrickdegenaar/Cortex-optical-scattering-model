% =========================================================================
%  EmissionPhaseArray.m  —  LED initial emission (theta, phi) LUT
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Returns a 2 x N look-up table of initial photon emission angles for a
%  surface-emitting LED of selectable angular profile.
%      Row 1:  polar angle theta in [-pi/2, pi/2]  (relative to surface
%              normal; negative side mirrors positive)
%      Row 2:  azimuthal angle phi  in [0, 2*pi)
%
%  PER-THETA SAMPLING WEIGHTS:
%    isotropic         w(theta) ~ |sin theta|
%    lambertian        w(theta) ~ cos(theta) * |sin theta|        (see [1])
%    tyndall           empirical top + sidewall LED pattern
%                      * |sin theta| * |cos theta|
%    pseudo collimated w(theta) ~ cos^50(theta) * |sin theta|
%    collimated        delta(theta - 0)  (all photons along the normal)
%
%  The cos(theta) in the Lambertian case is the Lambertian intensity law
%  I(theta) = I_0 cos(theta) (Lambert's cosine law, [1,2]).  The sin(theta)
%  Jacobian comes from the solid-angle element d Omega = sin theta d theta
%  d phi (needed to turn a per-steradian intensity into a per-theta PDF).
%
%  References:
%    [1] Lambert, J.H. (1760), "Photometria, sive de mensura et gradibus
%        luminis, colorum et umbrae", Augsburg.   (Lambert's cosine law.)
%    [2] McCluney, W.R. (1994), "Introduction to Radiometry and Photometry",
%        Artech House, Sec. 2.2.
%    [3] Schubert, E.F. (2006), "Light-Emitting Diodes", 2nd ed., Cambridge
%        University Press, Ch. 5  (LED radiation patterns).
% =========================================================================

function result = EmissionPhaseArray(SourceType, Angle_step)

    Arc_step  = Angle_step * pi/180;
    num_angle = ceil(pi / Arc_step);
    theta0    = -pi/2 : pi/num_angle : pi/2;

    switch SourceType
        case 'isotropic'
            % Uniform on the upward hemisphere
            prob = abs(sin(theta0));

        case 'lambertian'
            % I(theta) = I_0 cos(theta)  -->  PDF in theta has extra sin(theta)
            prob = cos(theta0) .* abs(sin(theta0));

        case 'tyndall'
            % Empirical surface-emitting LED radiation pattern (top lobe
            % + sidewall lobes).  Fitted, not derived from first principles.
            angles1  = -15:Angle_step:15;
            piAngles = angles1 * pi/180;
            topEm    = cos(abs(piAngles)/2) - 1.05*cos(pi/8);
            topEm    = topEm / max(topEm);

            angles2  = (15 + Angle_step):Angle_step:90;
            sideEm   = 1 ./ (angles2.^0.25);
            sideEm   = sideEm - min(sideEm);
            sideEm   = sideEm * cos(pi/4) / max(sideEm);

            radDist_Tyn = [fliplr(sideEm), topEm, sideEm];
            if numel(radDist_Tyn) ~= numel(theta0)
                radDist_Tyn = interp1( ...
                    linspace(-pi/2, pi/2, numel(radDist_Tyn)), ...
                    radDist_Tyn, theta0, 'linear', 0);
            end
            prob = radDist_Tyn .* abs(sin(theta0)) .* abs(cos(theta0));

        case 'pseudo collimated'
            % Tight forward cone modelled as cos^N theta with N=50
            prob = (cos(theta0).^50) .* abs(sin(theta0));

        case 'collimated'
            % Delta at theta=0: all photons launched along the surface normal
            num_phai    = ceil(2*pi/Arc_step);
            theta_array = zeros(1, num_phai+1);
            phai_array  = zeros(1, num_phai+1);
            result      = [theta_array; phai_array];
            return

        otherwise
            error('EmissionPhaseArray:Unknown', ...
                  'Unknown SourceType "%s"', SourceType);
    end

    % --- LUT expansion (cumsum, no while-loop stalling on equal bins) ---
    prob = prob ./ sum(prob);
    prob = ceil(50000 .* prob);

    Oprob  = cumsum(prob);
    totalN = Oprob(end);

    theta_array = zeros(1, totalN);
    theta_array(1:Oprob(1)) = theta0(1);
    for u = 2:length(theta0)
        if prob(u) > 0
            theta_array(Oprob(u-1)+1:Oprob(u)) = theta0(u);
        end
    end

    % Azimuth phi: uniform on [0, 2*pi)
    phai_array = (1:totalN) * 2*pi / totalN;

    result = [theta_array; phai_array];
end
