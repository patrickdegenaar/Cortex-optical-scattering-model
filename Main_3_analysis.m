% =========================================================================
%  Main_3_analysis.m  —  Post-process MC fluence map into irradiance
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Loads the 3D photon-deposition cube produced by Main_2 and
%      1. Smooths it with a small 3D Gaussian kernel (suppresses voxel noise).
%      2. Identifies the source plane by integrating each XZ/YZ row sum and
%         taking the peak (typically slice idx == 2 due to convolution edge
%         effects at z = 1).
%      3. Normalises so that the integrated weight at the source plane
%         equals 1 -> per-voxel value is a fraction of input power.
%      4. Divides by voxel cross-section (mm^2) -> per-voxel irradiance.
%      5. Renders linear and log10 XY / XZ / YZ slices.
%      6. Computes per-z-layer total power as an energy-conservation check.
%
%  CALIBRATION  (mW/mm^2):
%    The output  calibratedImage_xyz  is a *fraction of source power per
%    mm^2*.  Multiply by the source power (mW) at the emitting face to
%    obtain absolute irradiance.  For optogenetics, compare to the opsin
%    activation threshold (typically ~1 mW/mm^2 at 470 nm for ChR2; lower
%    for red-shifted opsins such as Chrimson at 617 nm) [3].
%
%  References:
%    [1] Wang, L., Jacques, S.L., Zheng, L. (1995), "MCML—Monte Carlo
%        modeling of light transport in multi-layered tissues", Computer
%        Methods and Programs in Biomedicine 47, 131-146.
%    [2] Jacques, S.L. & Pogue, B.W. (2008), "Tutorial on diffuse light
%        transport", J. Biomed. Opt. 13, 041302.
%    [3] Klapoetke, N.C. et al. (2014), "Independent optical excitation of
%        distinct neural populations", Nature Methods 11, 338-346.
%        (Chrimson red-shifted opsin, activation thresholds.)
% =========================================================================

close all
clear;
clc;

% Mac / POSIX:
load('scatterData/lambertian_D200um.mat');
% Windows: load('scatterData\\lambertian_D200um.mat');

%% simulation parameters (must match Main_2)
xMin = -750;  xMax = 750;
yMin = -750;  yMax = 750;
zMin =    0;  zMax = 1500;
Z_Offset = 100;

intensityThreshold = 1e-5;
intensityScale     = 1/intensityThreshold;

%--------------------------------------------------------------------------
%% Box dimensions
%--------------------------------------------------------------------------
[simBoxSizeX, simBoxSizeY, simBoxSizeZ] = size(scatterImage_xyz);

x_div = (xMax - xMin)/simBoxSizeX;       % um per voxel
y_div = (yMax - yMin)/simBoxSizeY;
z_div = (zMax - zMin)/simBoxSizeZ;

%--------------------------------------------------------------------------
%% 3D Gaussian smoothing
%
%  Convolution with a small isotropic Gaussian kernel removes the
%  pixel-by-pixel Monte-Carlo noise without distorting the underlying
%  fluence distribution provided the kernel width <= local correlation
%  length of the noise.  (See e.g. [2] for noise control in tissue MC.)
%--------------------------------------------------------------------------
gaussianSize = 11;
GaussianDist = 0.75;

kernel3           = gaussianGenerator3D(gaussianSize, GaussianDist);
smoothedImage_xyz = convn(scatterImage_xyz, kernel3, 'same');

%--------------------------------------------------------------------------
%% Normalisation step 1: peak normalisation (for display)
%--------------------------------------------------------------------------
smoothedImage_xyz = smoothedImage_xyz / max(smoothedImage_xyz(:));

%--------------------------------------------------------------------------
%% Plane extraction
%--------------------------------------------------------------------------
centrepointX = round(simBoxSizeX/2);
centrepointY = round(simBoxSizeY/2);

XYplane_z = round(simBoxSizeZ/2);
XZplane_y = round(simBoxSizeY/2);
YZplane_x = round(simBoxSizeX/2);

smoothedImage_yz = squeeze(smoothedImage_xyz(YZplane_x, :, :));
smoothedImage_xz = squeeze(smoothedImage_xyz(:, XZplane_y, :));
smoothedImage_xy = smoothedImage_xyz(:, :, XYplane_z);

%--------------------------------------------------------------------------
%% Identify the source plane (max column-sum)
%
%  Edge effects of the smoothing convolution can push the apparent
%  brightest z-slice off the surface by 1 voxel; we pick whichever slice
%  carries the most integrated weight after smoothing.
%--------------------------------------------------------------------------
sumVals = zeros(1, simBoxSizeY);
for i = 1:simBoxSizeY
    sumVals(i) = sum(smoothedImage_yz(:, i));
end
[~, idx] = max(sumVals);

%--------------------------------------------------------------------------
%% Normalisation step 2:  fraction-of-source-power per voxel
%
%  Integrated counts at the source plane == total input weight (~num_photon).
%  Dividing the whole cube by this number turns each voxel into a fraction
%  of the source power.  Dividing again by voxel area (mm^2) gives
%  irradiance (mW/mm^2 / mW_source).
%--------------------------------------------------------------------------
maxPlane = smoothedImage_xyz(:, :, idx);
sumPlane = sum(maxPlane(:));

calibratedImage_xyz = smoothedImage_xyz / sumPlane;

voxel_area = (x_div/1000) * (y_div/1000);                  % mm^2
calibratedImage_xyz = calibratedImage_xyz / voxel_area;     % per mm^2

%--------------------------------------------------------------------------
%% Re-extract planes from the calibrated cube
%--------------------------------------------------------------------------
calibratedImage_yz = squeeze(calibratedImage_xyz(YZplane_x, :, :));
calibratedImage_xz = squeeze(calibratedImage_xyz(:, XZplane_y, :));
calibratedImage_xy = calibratedImage_xyz(:, :, XYplane_z);

logNormalisedImage_yz = log10(calibratedImage_yz);
logNormalisedImage_xz = log10(calibratedImage_xz);
logNormalisedImage_xy = log10(calibratedImage_xy);

%--------------------------------------------------------------------------
%% Display axes
%--------------------------------------------------------------------------
xAxis = xMin:x_div:xMax - x_div;
yAxis = yMin:y_div:yMax - y_div;
zAxis = zMin:z_div:zMax - z_div;

%--------------------------------------------------------------------------
%% Linear-scale plots
%--------------------------------------------------------------------------
figure; set(gcf, 'Color', 'white');

subplot(2,2,1);
pcolor(zAxis, yAxis, calibratedImage_yz);
shading interp; axis equal tight; title('YZ plane');

subplot(2,2,2);
pcolor(zAxis, xAxis, calibratedImage_xz);
shading interp; axis equal tight; title('XZ plane');

subplot(2,2,3);
pcolor(xAxis, yAxis, calibratedImage_xy);
shading interp; axis equal tight; title('XY plane');

subplot(2,2,4);
axis off; colorbar; title('Linear irradiance (mW/mm^2 per mW_{source})');

%--------------------------------------------------------------------------
%% Log-scale plots (3 decades)
%--------------------------------------------------------------------------
figure; set(gcf, 'Color', 'white');

subplot(2,2,1);
pcolor(zAxis, yAxis, logNormalisedImage_yz);
shading interp; axis equal tight; title('YZ plane'); clim([-3 0]);

subplot(2,2,2);
pcolor(zAxis, xAxis, logNormalisedImage_xz);
shading interp; axis equal tight; title('XZ plane'); clim([-3 0]);

subplot(2,2,3);
pcolor(xAxis, yAxis, logNormalisedImage_xy);
shading interp; axis equal tight; title('XY plane'); clim([-3 0]);

subplot(2,2,4);
axis off; colorbar; clim([-3 0]); title('log_{10} irradiance (3 decades)');

%--------------------------------------------------------------------------
%% Save calibrated cube
%--------------------------------------------------------------------------
if ~exist('GeneratedData', 'dir'); mkdir('GeneratedData'); end
save(fullfile('GeneratedData','lambertian_D200um_calibrated.mat'), ...
     'calibratedImage_xyz', '-v7.3');

%--------------------------------------------------------------------------
%% Write rescaled 8-bit JPGs (1 pixel = 1 um)
%--------------------------------------------------------------------------
[yDim, zDim] = size(calibratedImage_yz);
[xDim, yDim] = size(calibratedImage_xy);

xDim_px = xDim * x_div;
yDim_px = yDim * y_div;
zDim_px = zDim * z_div;

writeIm_xy = imresize(calibratedImage_xy, [xDim_px, yDim_px]);
writeIm_xz = imresize(calibratedImage_xz, [xDim_px, zDim_px]);
writeIm_yz = imresize(calibratedImage_yz, [yDim_px, zDim_px]);

writeIm_Log_xy = imresize(logNormalisedImage_xy, [xDim_px, yDim_px]);
writeIm_Log_xz = imresize(logNormalisedImage_xz, [xDim_px, zDim_px]);
writeIm_Log_yz = imresize(logNormalisedImage_yz, [yDim_px, zDim_px]);

writeIm_xy = uint8(writeIm_xy * 255/max(writeIm_xy(:)));
writeIm_xz = uint8(writeIm_xz * 255/max(writeIm_xz(:)));
writeIm_yz = uint8(writeIm_yz * 255/max(writeIm_yz(:)));

writeIm_Log_xy = writeIm_Log_xy - min(writeIm_Log_xy(:));
writeIm_Log_xz = writeIm_Log_xz - min(writeIm_Log_xz(:));
writeIm_Log_yz = writeIm_Log_yz - min(writeIm_Log_yz(:));
writeIm_Log_xy = uint8(writeIm_Log_xy * 255/max(writeIm_Log_xy(:)));
writeIm_Log_xz = uint8(writeIm_Log_xz * 255/max(writeIm_Log_xz(:)));
writeIm_Log_yz = uint8(writeIm_Log_yz * 255/max(writeIm_Log_yz(:)));

imwrite(writeIm_xy, 'linearXY.jpg');
imwrite(writeIm_xz, 'linearXZ.jpg');
imwrite(writeIm_yz, 'linearYZ.jpg');
imwrite(writeIm_Log_xy, 'logXY.jpg');
imwrite(writeIm_Log_xz, 'logXZ.jpg');
imwrite(writeIm_Log_yz, 'logYZ.jpg');

%--------------------------------------------------------------------------
%% Per-layer total power (energy-conservation check)
%
%  In a purely absorbing medium, integrating over each constant-z slab
%  gives the total power crossing that depth.  For tissue with scattering,
%  this is the *fluence-area integral*, which should be monotonically
%  non-increasing with z because of absorption + backscatter losses.
%--------------------------------------------------------------------------
layers = zeros(1, simBoxSizeZ);
for layer = 1:simBoxSizeZ
    layers(layer) = sum(sum(calibratedImage_xyz(:, :, layer))) * voxel_area;
end
figure; plot((1:simBoxSizeZ)*z_div, layers);
xlabel('z (\mum)'); ylabel('Layer-integrated power / source power');
title('Energy-conservation check'); grid on;
