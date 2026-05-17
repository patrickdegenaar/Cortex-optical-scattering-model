% =========================================================================
%  Main_2_scatterSimulation.m  —  3D Monte-Carlo photon-transport simulation
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Standard tissue-optics Monte-Carlo following the framework of
%  Wang/Jacques MCML [1] and Prahl et al. [2]:
%
%    For each photon:
%      1. Pick a source position on the LED face (Source_Table).
%      2. Pick an initial emission direction (Emitter_angles).
%      3. Loop:
%           a. Sample free path s from p(s) = (1/mfp) exp(-s/mfp)   [3]
%           b. Move photon by (s, theta, phi).
%           c. Apply continuous absorption to weight:
%                W <- W * exp(-mu_a * s)                            [1]
%           d. Choose scattering kernel (Mie / Rayleigh / isotropic),
%              draw local deflection (delta_theta, delta_phi).
%           e. Rotate local direction back to lab frame using R^T (R the
%              proper rotation aligning lab z with current direction).
%           f. Deposit W into the voxel at the new position.
%      4. Terminate on box exit or W < intensityThreshold.
%
%  ROTATION MATRIX  (proper, right-handed, det = +1):
%      R = [  cos(theta)*cos(phi)  cos(theta)*sin(phi)  -sin(theta);
%             -sin(phi)            cos(phi)               0        ;
%             sin(theta)*cos(phi)  sin(theta)*sin(phi)   cos(theta)]
%  rows = local basis (theta-hat, phi-hat, r-hat) in the lab frame.
%  d_lab = R' * d_local  because R is orthogonal (R^-1 = R').
%  Polar-angle recovery uses acos(d_z) and atan2(d_y, d_x) — non-singular
%  everywhere (cf. the standard MCML direction-cosine update which has a
%  removable singularity at |u_z| -> 1).
%
%  PHASE-FUNCTION SAMPLING:  delta_theta is drawn from the LUT built by
%  Mie_scat_angle_array / Rayleigh_scat_angle_array / Even_scat_angle_array.
%  Each LUT already includes the sin(theta) solid-angle Jacobian, so
%  rejection-free uniform-index sampling gives the correct distribution.
%
%  ABSORPTION SCHEME:  continuous "weight" accounting (Wang/Jacques) rather
%  than the binary albedo-rejection variant; statistically equivalent in
%  expectation but lower variance per photon.
%
%  References:
%    [1] Wang, L., Jacques, S.L., Zheng, L. (1995), "MCML—Monte Carlo
%        modeling of light transport in multi-layered tissues", Computer
%        Methods and Programs in Biomedicine 47, 131-146.
%    [2] Prahl, S.A., Keijzer, M., Jacques, S.L., Welch, A.J. (1989),
%        "A Monte Carlo Model of Light Propagation in Tissue", SPIE
%        Institute Series IS 5, 102-111.
%    [3] Ishimaru, A. (1978), "Wave Propagation and Scattering in Random
%        Media", Academic Press, Ch. 7  (Beer-Lambert / free path PDF).
%    [4] Jacques, S.L. & Pogue, B.W. (2008), "Tutorial on diffuse light
%        transport", J. Biomed. Opt. 13, 041302.
% =========================================================================

clear all
clc

n = 1;
dataFile = struct('fileName', 'lambertian_D200um');
% dataFile = struct('fileName', 'tyndall');
% dataFile = struct('fileName', 'pseudo collimated');
% dataFile = struct('fileName', 'collimated');
% dataFile = struct('fileName', 'isotropic');

% On a Windows system, comment-out the POSIX line below and use:
%   file = strcat('settingsData\Settings_', dataFile.fileName);
file = strcat('settingsData/Settings_lambertian_D200um.mat');

load(file);

% -------------------------------------------------------------------------
%  Scattering-type mixture
% -------------------------------------------------------------------------
%   Prob_scattering  — probability a scattering event occurs at all
%                      (1.0 = always scatter; 0.0 = pure absorption)
%   Coeff_*          — fractional contribution of each kernel
% -------------------------------------------------------------------------
Prob_scattering = 1.0;
Coeff_Mie = 0.95;                       % brain at red wavelengths is Mie-dominated
Coeff_Ryl = 0.05;                       % small Rayleigh tail (small organelles)
Coeff_Ist = 1 - Coeff_Mie - Coeff_Ryl;

% Decision thresholds for kernel selection on a single uniform draw:
Criterion1 = 1 - Prob_scattering;                       % no scattering
Criterion2 = Criterion1 + Prob_scattering * Coeff_Mie;  % Mie
Criterion3 = Criterion2 + Prob_scattering * Coeff_Ryl;  % Rayleigh
Criterion4 = 1;                                          % isotropic

% -------------------------------------------------------------------------
%  Simulation parameters
% -------------------------------------------------------------------------
num_photon         = 5e5;
intensityThreshold = 1e-5;
intensityScale     = 1/intensityThreshold;

% Simulation box (um)
xMin = -1500;   xMax = 1500;
yMin = -1500;   yMax = 1500;
zMin =     0;   zMax = 1500;
Z_Offset = 100;

% Image-processing knobs (used downstream by Main_3)
gaussianSize = 11;        % convolution kernel size (odd)
GaussianDist = 0.75;      % stddev of Gaussian kernel

% Voxel grid
simBoxSizeX = 600;
simBoxSizeY = 600;
simBoxSizeZ = 600;

x_div = (xMax - xMin)/simBoxSizeX;       % um per voxel
y_div = (yMax - yMin)/simBoxSizeY;
z_div = (zMax - zMin)/simBoxSizeZ;

xAxis = xMin:x_div:xMax - x_div;
yAxis = yMin:y_div:yMax - y_div;
zAxis = zMin:z_div:zMax - z_div;

% Output 3D fluence map — pre-filled with a small "floor" value so log-plot
% sweeps over the whole volume even before photons arrive.
scatterImage_xyz(1:simBoxSizeX, 1:simBoxSizeY, 1:simBoxSizeZ) = 1e-5;

% Mean free path used inside the loop (LUT-sampled per step)
dt_r = mean(Photon_paths);

% =========================================================================
%  Monte-Carlo loop
% =========================================================================
tic
for i = 1:num_photon

    % ---- 1. Pick source location (uniformly from Source_Table) ----------
    SourceOrder = ceil(size_SrcTable(1,1) * rand);
    xi = Source_Table(SourceOrder, 1);
    yi = Source_Table(SourceOrder, 2);
    zi = Source_Table(SourceOrder, 3);

    x(1,1) = xi;  y(1,1) = yi;  z(1,1) = zi;
    j      = 1;

    % ---- 2. Pick initial direction from emission LUT --------------------
    theta(1,1) = Emitter_angles_Lat(ceil(size_inital_array(1,2) * rand));
    phai (1,1) = Emitter_angles_Azu(ceil(size_inital_array(1,2) * rand));

    % Fold theta in (-pi/2, 0) into (0, pi/2) with phi shifted by pi
    if theta(1,1) < 0
        theta(1,1) = -theta(1,1);
        phai (1,1) =  phai(1,1) + pi;
    end

    % ---- 3. Propagation loop --------------------------------------------
    intensity = 1;
    while intensity >= intensityThreshold

        j = j + 1;

        % --- (a) Sample free path (Beer-Lambert) -------------------------
        det_r2 = Photon_paths(randi(numel(Photon_paths)));

        % --- (b) Move photon ---------------------------------------------
        z(1,j) = z(1,j-1) + det_r2 * cos(theta(1,j-1));
        x(1,j) = x(1,j-1) + det_r2 * sin(theta(1,j-1)) * cos(phai(1,j-1));
        y(1,j) = y(1,j-1) + det_r2 * sin(theta(1,j-1)) * sin(phai(1,j-1));

        % --- Boundary conditions: discard photons leaving the box --------
        if z(1,j) <= zMin || z(1,j) >= zMax || x(1,j) <= xMin ...
                || x(1,j) >= xMax || y(1,j) <= yMin || y(1,j) >= yMax
            intensity = 0;
            break
        end

        % --- (c) Continuous absorption (Beer-Lambert over free path) -----
        %   I -> I * exp(-mu_a * s)
        %   Applies to every path segment, scattered or not.
        ls_factor = exp(-coeff_absorb * det_r2);
        intensity = intensity * ls_factor;

        % --- (d) Pick scattering kernel ----------------------------------
        scattering_switch = rand;

        if scattering_switch <= Criterion1
            % No scattering this step
            theta(1,j) = theta(1,j-1);
            phai (1,j) = phai (1,j-1);
        else
            if scattering_switch <= Criterion2
                delta_theta = scat_angles_lateral_Mie(ceil(size_angle_array_Mie(1,2)*rand));
                delta_phai  = scat_angles_azimuth_Mie(ceil(size_angle_array_Mie(1,2)*rand));
            elseif scattering_switch <= Criterion3
                delta_theta = scat_angles_lateral_Ryl(ceil(size_angle_array_Ryl(1,2)*rand));
                delta_phai  = scat_angles_azimuth_Ryl(ceil(size_angle_array_Ryl(1,2)*rand));
            else
                delta_theta = scat_angles_lateral_Ist(ceil(size_angle_array_Ist(1,2)*rand));
                delta_phai  = scat_angles_azimuth_Ist(ceil(size_angle_array_Ist(1,2)*rand));
            end

            % --- (e) Rotation: local -> lab frame ------------------------
            % Build the proper rotation matrix R that maps the lab frame
            % to the local frame in which the photon travels along z_loc.
            % Rows of R are (theta-hat, phi-hat, r-hat) of the spherical
            % triad evaluated at the current (theta, phi).  See header
            % for the algebra.
            if delta_theta == 0
                theta(1,j) = theta(1,j-1);
                phai (1,j) = phai (1,j-1);
            else
                ct = cos(theta(1,j-1));  st = sin(theta(1,j-1));
                cp = cos(phai (1,j-1));  sp = sin(phai (1,j-1));

                R = [ ct*cp,  ct*sp, -st ;
                      -sp,    cp,     0  ;
                      st*cp,  st*sp,  ct ];

                % Direction in local frame after scattering by (dtheta,dphi)
                dv_loc = [ sin(delta_theta)*cos(delta_phai) ;
                           sin(delta_theta)*sin(delta_phai) ;
                           cos(delta_theta)                 ];

                % Local -> lab.  R orthogonal => inverse is transpose.
                dv_lab = R.' * dv_loc;
                dv_lab = dv_lab / norm(dv_lab);            % round-off guard

                % Recover (theta, phi) from direction cosines.
                % atan2 is well-defined everywhere — replaces the older
                % cos_phi = dv(1)/sqrt(1-dv(3)^2) formulation that blows
                % up near the poles.
                theta(1,j) = acos( max(-1, min(1, dv_lab(3))) );
                phai (1,j) = mod(  atan2(dv_lab(2), dv_lab(1)),  2*pi );
            end
        end

        % --- (f) Deposit weight into voxel grid --------------------------
        rescaledZ = ceil( (z(1,j) - zMin)/z_div );
        rescaledX = ceil( (xMax   - x(1,j))/x_div );
        rescaledY = ceil( (yMax   - y(1,j))/y_div );

        if rescaledX >= 1 && rescaledX <= simBoxSizeX && ...
           rescaledY >= 1 && rescaledY <= simBoxSizeY && ...
           rescaledZ >= 1 && rescaledZ <= simBoxSizeZ
            scatterImage_xyz(rescaledX, rescaledY, rescaledZ) = ...
                scatterImage_xyz(rescaledX, rescaledY, rescaledZ) + intensity;
        end
    end
end
simTime = toc;
disp(strcat('Monte Carlo loop time =  ', num2str(simTime), ' s'));

% =========================================================================
%  Prepare 2D slices for quick visualisation
% =========================================================================
centrepointX = round(simBoxSizeX/2);
centrepointY = round(simBoxSizeY/2);
centrepointZ = round(simBoxSizeZ/2);

scatterImage_xz(1:simBoxSizeX, 1:simBoxSizeZ) = scatterImage_xyz(:, centrepointY, :);
scatterImage_yz(1:simBoxSizeY, 1:simBoxSizeZ) = scatterImage_xyz(centrepointX, :, :);

% XY slice near the emitter face (z near Z_Offset)
z_idx = max(1, round((Z_Offset - zMin)/z_div));
scatterImage_xy(1:simBoxSizeX, 1:simBoxSizeY) = scatterImage_xyz(:, :, z_idx);

% =========================================================================
%  Plot
% =========================================================================
figure;
min_scalebar = -3;
maxVal       = num_photon / 10000 / pi;    % display normalisation

subplot(2,2,1);
pcolor(zAxis, yAxis, log10(scatterImage_yz./maxVal));
shading interp; axis equal tight; clim([min_scalebar 0]); title('YZ plane');

subplot(2,2,2);
pcolor(zAxis, xAxis, log10(scatterImage_xz./maxVal));
shading interp; axis equal tight; clim([min_scalebar 0]); title('XZ plane');

subplot(2,2,3);
pcolor(log10(scatterImage_xy./maxVal));
shading interp; axis equal tight; clim([min_scalebar 0]); title('XY plane');

subplot(2,2,4);
axis off; colorbar; clim([min_scalebar 0]); title(dataFile.fileName);

% =========================================================================
%  Save
% =========================================================================
file = strcat('ScatterData_', dataFile.fileName, '.mat');

% Check to see if a scatterData directory exists. If not insert. 
if ~exist('scatterData', 'dir'); mkdir('scatterData'); end

% Mac / POSIX
save(strcat('scatterData/', file), 'scatterImage_xyz', '-v7.3');
save('scatterData/lambertian_D200um.mat', 'scatterImage_xyz', '-v7.3');

% Windows alternatives:
% save(strcat('scatterData\', file), 'scatterImage_xyz', '-v7.3');
% save('scatterData\lambertian_D200um.mat', 'scatterImage_xyz', '-v7.3');
