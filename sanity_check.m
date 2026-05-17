% =========================================================================
%  sanity_check.m  —  Statistical sanity test of the LUT generators
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Runs Main_1_settingsGenerator, draws N samples from each look-up table
%  (free path, isotropic theta, Rayleigh theta, Mie theta, Lambertian
%  emission theta), and overlays the analytic / theoretical PDF for
%  visual verification.  Saves sanity_check.png.
%
%  Each LUT should match its analytic PDF within Monte-Carlo binning noise.
%  Discrepancies at low-density tails are expected and harmless (integer-
%  quantisation of the cumulative table).
%
%  References:
%    [1] Wang, L., Jacques, S.L., Zheng, L. (1995), MCML, CMPB 47, 131-146.
%    [2] Prahl, S.A. et al. (1989), SPIE Inst. Series IS 5.
%    [3] Bohren, C.F. & Huffman, D.R. (1983), "Absorption and Scattering of
%        Light by Small Particles", Wiley.
% =========================================================================

clear; clc; close all;
addpath(pwd);

% ---- 1. Run the settings generator --------------------------------------
Main_1_settingsGenerator

% ---- 2. Sample each LUT -------------------------------------------------
N = 200000;
rng(42);

Lp        = numel(Photon_paths);
s_samples = Photon_paths(randi(Lp, 1, N));

LI     = size(angle_array_Ist, 2);
th_iso = angle_array_Ist(1, randi(LI, 1, N));

LR     = size(angle_array_Ryl, 2);
th_ray = angle_array_Ryl(1, randi(LR, 1, N));

LM     = size(angle_array_Mie, 2);
th_mie = angle_array_Mie(1, randi(LM, 1, N));

LE      = size(Emitter_angles, 2);
th_emit = Emitter_angles(1, randi(LE, 1, N));

% ---- 3. Analytic / theoretical PDFs -------------------------------------
mfp     = 1/coeff_extinct;
s_grid  = linspace(0, 6*mfp, 400);
pdf_exp = (1/mfp) * exp(-s_grid/mfp);                       % Beer-Lambert [1]

th_grid = linspace(0, pi, 200);
pdf_iso = 0.5 * sin(th_grid);                                % uniform-on-sphere
pdf_ray = (3/8) * (1 + cos(th_grid).^2) .* sin(th_grid);     % Rayleigh
pdf_ray = pdf_ray / trapz(th_grid, pdf_ray);

nTH         = length(th_grid);
pdf_mie_raw = zeros(1, nTH);
for u = 1:nTH
    pdf_mie_raw(u) = Mie_scat_pattern(m, a, cos(th_grid(u)));
end
pdf_mie = pdf_mie_raw .* sin(th_grid);                      % apply sin Jacobian
pdf_mie = pdf_mie / trapz(th_grid, pdf_mie);

the_grid = linspace(-pi/2, pi/2, 200);
pdf_emit = cos(the_grid) .* abs(sin(the_grid));              % Lambertian
pdf_emit = pdf_emit / trapz(the_grid, pdf_emit);

% ---- 4. Plot ------------------------------------------------------------
figure('Position', [50 50 1200 800]);

% Octave/MATLAB-portable histogram helper
binPDF = @(x, nb, lo, hi) deal( ...
    lo + (hi-lo)*((0.5:nb-0.5)/nb), ...
    histc(x, lo + (hi-lo)*(0:nb)/nb)(1:nb) / numel(x) / ((hi-lo)/nb) );

subplot(2,3,1);
[cs, hs] = binPDF(s_samples, 60, 0, 6*mfp);
bar(cs, hs, 1.0, 'FaceColor', [0.7 0.7 0.9], 'EdgeColor', 'none'); hold on;
plot(s_grid, pdf_exp, 'r-', 'LineWidth', 2);
xlabel('free path s (\mum)'); ylabel('pdf');
title(sprintf('Free-path (mfp = %.2f \\mum)', mfp));
legend('LUT sample', 'exp(-s/mfp)/mfp', 'Location', 'northeast'); grid on;

subplot(2,3,2);
[cs, hs] = binPDF(th_iso, 60, 0, pi);
bar(cs, hs, 1.0, 'FaceColor', [0.7 0.9 0.7], 'EdgeColor', 'none'); hold on;
plot(th_grid, pdf_iso, 'r-', 'LineWidth', 2);
xlabel('\theta (rad)'); ylabel('pdf');
title('Isotropic   p(\theta)=½sin\theta'); xlim([0 pi]); grid on;

subplot(2,3,3);
[cs, hs] = binPDF(th_ray, 60, 0, pi);
bar(cs, hs, 1.0, 'FaceColor', [0.9 0.8 0.7], 'EdgeColor', 'none'); hold on;
plot(th_grid, pdf_ray, 'r-', 'LineWidth', 2);
xlabel('\theta (rad)'); ylabel('pdf');
title('Rayleigh   p(\theta)\propto(1+cos^2\theta)sin\theta'); xlim([0 pi]); grid on;

subplot(2,3,4);
[cs, hs] = binPDF(th_mie, 90, 0, pi);
bar(cs, hs, 1.0, 'FaceColor', [0.9 0.7 0.9], 'EdgeColor', 'none'); hold on;
plot(th_grid, pdf_mie, 'r-', 'LineWidth', 2);
xlabel('\theta (rad)'); ylabel('pdf');
title(sprintf('Mie  (x=%.2f, m=%.3f+%.1ei)', a, real(m), imag(m)));
xlim([0 pi]); grid on;

subplot(2,3,5);
[cs, hs] = binPDF(th_emit, 80, -pi/2, pi/2);
bar(cs, hs, 1.0, 'FaceColor', [0.7 0.9 0.9], 'EdgeColor', 'none'); hold on;
plot(the_grid, pdf_emit, 'r-', 'LineWidth', 2);
xlabel('\theta (rad)'); ylabel('pdf');
title('Lambertian emission cos\theta\cdot|sin\theta|');
xlim([-pi/2 pi/2]); grid on;

% Numerical summary panel
subplot(2,3,6); axis off;
mean_mie = trapz(th_grid, th_grid.*pdf_mie);
txt = {
    '\bf Sanity-check summary'
    sprintf('Free-path:   mean(LUT)=%.3f um   mean(th)=%.3f um', mean(s_samples), mfp)
    sprintf('Isotropic:   mean(LUT)=%.3f rad  mean(th)=%.3f rad', mean(th_iso), pi/2)
    sprintf('Rayleigh:    mean(LUT)=%.3f rad  mean(th)=%.3f rad', mean(th_ray), pi/2)
    sprintf('Mie:         mean(LUT)=%.3f rad  mean(th)=%.3f rad', mean(th_mie), mean_mie)
    ''
    sprintf('Bulk optical (current settings):')
    sprintf('   mu_ext = %.3f mm^{-1}', coeff_extinct*1e3)
    sprintf('   mu_abs = %.3f mm^{-1}', coeff_absorb*1e3)
    sprintf('   mu_sca = %.3f mm^{-1}', (coeff_extinct-coeff_absorb)*1e3)
    sprintf('   mfp    = %.3f um',      mfp)
};
text(0, 0.95, txt, 'Units', 'normalized', 'VerticalAlignment', 'top', ...
     'FontName', 'Helvetica', 'FontSize', 10);

saveas(gcf, 'sanity_check.png');
fprintf('\nSaved: sanity_check.png\n');
