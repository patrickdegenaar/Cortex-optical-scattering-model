% =========================================================================
%  SourceArray.m  —  Spatial sampling of the LED emitting surface
%
%  Reviewed by Claude Opus 4.7 for correctness on 16 May 2026.
%
%  Returns a probability-weighted look-up table of (x,y,z=0) point sources
%  distributed across the circular face of an LED of diameter d_src.  The
%  spatial grid uses uniform radial sampling with arc-length-adapted
%  azimuthal sampling to keep the on-disk density approximately uniform
%  (the number of points per ring grows as ~r, so the per-area count is
%  roughly constant).
%
%  Per-point probability weights:
%    isotropic / lambertian / tyndall / pseudo collimated:
%        uniform weight 1/N    (intensity is uniform across the LED face)
%    collimated:
%        Gaussian profile centred at (0,0):
%            P(r) ~ exp( -coeff * r^2 ),   coeff = e^0.3 / (d_src/2)^2
%        (gives a beam-like spatial distribution rather than a flat top)
%
%  Note:  the angular distribution is handled separately by
%  EmissionPhaseArray.m; this script only sets the spatial origin.
%
%  References:
%    [1] McCluney, W.R. (1994), "Introduction to Radiometry and Photometry",
%        Artech House.
%    [2] Schubert, E.F. (2006), "Light-Emitting Diodes", 2nd ed., Cambridge
%        University Press.
% =========================================================================

function result = SourceArray(SourceType, d_src, r_step)

     Arc_step = r_step .* pi/2;          % azimuthal arc-length spacing seed

     k = 2;
     PSource(1,1) = 0;
     PSource(1,2) = 0;
     PSource(1,3) = 0;

     for r = r_step:r_step:d_src/2
         ang_step = Arc_step ./ r;
         for th = 0:ang_step:(2*pi - ang_step)
             PSource(k,1) = r .* cos(th);
             PSource(k,2) = r .* sin(th);
             PSource(k,3) = 0;
             k = k + 1;
         end
     end

     sz_psrc  = size(PSource);
     num_psrc = sz_psrc(1);

     %---- Per-point weights ------------------------------------------------
     switch SourceType
        case {'isotropic','lambertian','tyndall','pseudo collimated'}
            % flat-top emitter
            PSource(:,4) = 1/num_psrc;
        case 'collimated'
            % Gaussian spatial profile (beam-like)
            coeff = ((exp(1))^0.3) / ((d_src/2)^2);
            for i = 1:num_psrc
                PSource(i,4) = exp(-coeff*(PSource(i,1)^2 + PSource(i,2)^2));
            end
            PSource(:,4) = PSource(:,4) ./ sum(PSource(:,4));
     end

     %---- Expand to integer-occupancy LUT ---------------------------------
     prob = ceil( PSource(:,4) ./ min(PSource(:,4)) );
     Oprob(1) = prob(1);
     for u = 2:num_psrc
         Oprob(u) = Oprob(u-1) + prob(u);
     end

     u = 1;
     for v = 1:sum(prob)
         if v <= Oprob(1)
             source_array(v,1) = PSource(1,1);
             source_array(v,2) = PSource(1,2);
             source_array(v,3) = PSource(1,3);
         elseif v > Oprob(u) && v <= Oprob(u+1)
             source_array(v,1) = PSource(u+1,1);
             source_array(v,2) = PSource(u+1,2);
             source_array(v,3) = PSource(u+1,3);
             if v == Oprob(u+1)
                 u = u + 1;
             end
         end
     end

     result = source_array;
end
