% =========================================================================
%  Main_1_settingsGenerator.m  —  Pre-compute LUTs for the photon-transport
%                                  Monte-Carlo simulation
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Pipeline:
%      Main_1_settingsGenerator  -->  Main_2_scatterSimulation  -->  Main_3_analysis
%
%  Generates and saves to disk:
%      Source_Table              point-source positions (x,y,z) on LED face
%      Emitter_angles_Lat/Azu    initial (theta, phi) emission angles
%      Photon_paths              free-path-length LUT (exponential)
%      coeff_absorb, coeff_extinct   bulk absorption / extinction (um^-1)
%      scat_angles_*_{Ist,Ryl,Mie}   scattering-angle LUTs
%
%  PHYSICAL MODEL:
%    Bulk optical coefficients from single-particle Mie theory:
%        mu_ext  =  rho * pi (d/2)^2 * Q_ext(m,x)               (um^-1)
%        mu_abs  =  rho * pi (d/2)^2 * Q_abs(m,x)
%        mu_sca  =  mu_ext - mu_abs
%    with rho the number density of scatterers (particles um^-3), d the
%    scatterer diameter, x = pi d / (lambda/n_e) the size parameter, and
%    m = n_p/n_e the relative refractive index (complex).  The imaginary
%    part of n_p carries the absorption.  See Bohren & Huffman (1983).
%
%    For brain tissue at red wavelengths, target (Yaroslavsky 2002,
%    Jacques 2013):
%        617 nm  grey matter:  mu_a ~ 0.05-0.10 mm^-1,
%                              mu_s ~ 9-11 mm^-1,  g ~ 0.85-0.95
%
%  References:
%    [1] Bohren, C.F. & Huffman, D.R. (1983), "Absorption and Scattering of
%        Light by Small Particles", Wiley.
%    [2] Wang, L., Jacques, S.L., Zheng, L. (1995), "MCML—Monte Carlo
%        modeling of light transport in multi-layered tissues", Computer
%        Methods and Programs in Biomedicine 47, 131-146.
%    [3] Jacques, S.L. (2013), "Optical properties of biological tissues:
%        a review", Phys. Med. Biol. 58, R37-R61.
%    [4] Yaroslavsky, A.N. et al. (2002), "Optical properties of selected
%        native and coagulated human brain tissues in vitro in the visible
%        and near infrared spectral range", Phys. Med. Biol. 47, 2059-2073.
% =========================================================================

clear
clc
tic

% -------------------------------------------------------------------------
%  Settings
% -------------------------------------------------------------------------
SourceType = 'lambertian';     % 'isotropic'|'lambertian'|'tyndall'|
                               % 'pseudo collimated'|'collimated'

%  Part I  Light source spatial extent
diameter_source = 200;          % source diameter (um)
Radius_step     = 1;            % radial sampling (um)

% --- Wavelength-dependent optical parameters -----------------------------
%
%  Red light  (617 nm) — current selection
lamda = 0.617;                  % operating wavelength (um)
ne    = 1.355;                  % refractive index of medium (tissue)
np    = 1.47 + 0.10e-3*1i;      % refractive index of scatterer
                                % (Im part carries tissue absorption)

%  Blue light (470 nm) — alternative; uncomment to switch
% lamda = 0.47;
% ne    = 1.36;
% np    = 1.48 + 0.37e-3*1i;

% --- Scatterer population --------------------------------------------------
d     = 0.6;                    % scatterer diameter (um) — single-population
                                % approximation; real tissue has a distribution
rou_M = 2.5e-10;                % number concentration (mol/L); tuned to give
                                % mu_s in the brain literature range at 617 nm
rou   = rou_M * 6.02e8;         % --> particles / um^3   (= N_A x mol/L x 1e-15)

% --- Simulator sampling resolution ---------------------------------------
Angle_step = 0.5;               % angular bin width (deg)
path_step  = 0.5;               % free-path bin width (um)


% -------------------------------------------------------------------------
%  Part II  Mie efficiencies and bulk coefficients
% -------------------------------------------------------------------------
[m, a, pcs, pV] = Parameters(lamda, d, ne, np);  % m, x, pi*d^2/4, (4/3)pi*r^3

result_Mie = Mie(m, a);                          % [m' m" x Qext Qsca Qabs Qb g Qb/Qsca]
Q_ext = result_Mie(4);
Q_abs = result_Mie(6);

crossSection_ext = Q_ext * pcs;                  % um^2 per particle
crossSection_abs = Q_abs * pcs;

coeff_extinct = rou * crossSection_ext;          % mu_ext (um^-1)
coeff_absorb  = rou * crossSection_abs;          % mu_abs (um^-1)
coeff_scatter = coeff_extinct - coeff_absorb;    % mu_sca

% -------------------------------------------------------------------------
%  Part III  Photon-transport LUTs
% -------------------------------------------------------------------------
dt_r = 1 / coeff_extinct;       % mean free path (um) = 1/mu_t

% Scattering-angle LUTs (1-degree resolution, with sin(theta) Jacobian)
angle_array_Ist = Even_scat_angle_array(1);
angle_array_Ryl = Rayleigh_scat_angle_array(1);
angle_array_Mie = Mie_scat_angle_array(m, a);

scat_angles_lateral_Ist = angle_array_Ist(1,:);
scat_angles_azimuth_Ist = angle_array_Ist(2,:);
scat_angles_lateral_Ryl = angle_array_Ryl(1,:);
scat_angles_azimuth_Ryl = angle_array_Ryl(2,:);
scat_angles_lateral_Mie = angle_array_Mie(1,:);
scat_angles_azimuth_Mie = angle_array_Mie(2,:);

% Initial-phase LUT  (LED emission profile)
Emitter_angles     = EmissionPhaseArray(SourceType, Angle_step);
Emitter_angles_Lat = Emitter_angles(1,:);
Emitter_angles_Azu = Emitter_angles(2,:);

% Free-path LUT  (Beer-Lambert exponential)
Photon_paths = PathArray(path_step, dt_r);

% Source-position LUT
Source_Table = SourceArray(SourceType, diameter_source, Radius_step);

% Sizes (used by the MC loop in Main_2)
size_angle_array_Mie = size(angle_array_Mie);
size_angle_array_Ist = size(angle_array_Ist);
size_angle_array_Ryl = size(angle_array_Ryl);
size_inital_array    = size(Emitter_angles);
size_paths           = size(Photon_paths);
size_SrcTable        = size(Source_Table);

% --- Report ---------------------------------------------------------------
fprintf('lambda          = %.3f um\n', lamda);
fprintf('Q_ext, Q_abs    = %.4f, %.4e\n', Q_ext, Q_abs);
fprintf('mu_ext (um^-1)  = %.4e   (= %.3f mm^-1)\n', coeff_extinct, coeff_extinct*1e3);
fprintf('mu_abs (um^-1)  = %.4e   (= %.3f mm^-1)\n', coeff_absorb,  coeff_absorb*1e3);
fprintf('mu_sca (um^-1)  = %.4e   (= %.3f mm^-1)\n', coeff_scatter, coeff_scatter*1e3);
fprintf('g (anisotropy)  = %.4f\n', result_Mie(8));
fprintf('mfp    (um)     = %.3f\n', dt_r);

toc

% -------------------------------------------------------------------------
%  Save
% -------------------------------------------------------------------------
file = strcat('Settings_', SourceType, '_D', num2str(diameter_source), 'um.mat');

if ~exist('settingsData', 'dir'); mkdir('settingsData'); end
save(fullfile('settingsData', file), ...
     'Source_Table', 'Emitter_angles_Lat', 'Emitter_angles_Azu', ...
     'Photon_paths', 'coeff_absorb', 'coeff_extinct', ...
     'scat_angles_lateral_Mie', 'scat_angles_azimuth_Mie', ...
     'scat_angles_lateral_Ryl', 'scat_angles_azimuth_Ryl', ...
     'scat_angles_lateral_Ist', 'scat_angles_azimuth_Ist', ...
     'size_inital_array', 'size_SrcTable', 'size_paths', ...
     'size_angle_array_Mie', 'size_angle_array_Ryl', 'size_angle_array_Ist');
