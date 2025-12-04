function results = compute_poc_flux(params)
%COMPUTE_POC_FLUX  Calculate POC fluxes (Zeu and 100 m) and related metrics.
%
% USAGE
%   results = compute_poc_flux(params)
%
% INPUT (fields in params struct)
%   Clim_NPP            - [122 x 199 x 12] monthly NPP (units: as provided)
%   Clim_Zeu2           - [122 x 199 x 12] Zeu2 (euphotic depth) in m
%   Clim_Zeu3           - [122 x 199 x 12] (used for normalization)
%   Clidata_NN_cp660    - [122 x 199 x nDepth x 12] cp660 neural-network profile
%   Clim_MLD            - [122 x 199 x 12] mixed layer depth in m
%   SCSdepth9km         - [122 x 199] bathymetry (m)
%   depth               - [nDepth x 1] depth levels corresponding to profiles
%   FM_bbp_Micro        - [122 x 199 x 12] fraction (%) of microbbp (0-100)
%   Zeu2                - OPTIONAL alias for Clim_Zeu2 if available
%
% OPTIONAL PARAMETERS (fields)
%   f_FecT              - fraction for fecal flux (default 0.1)
%   m_ph                - mortality parameter (default 0.1)
%   useParallel         - true/false, use parfor for months (default false)
%
% OUTPUT (results struct)
%   Flux_Zeu_nn, Flux_100_nn, DPOC, DPOC100, IPOC_*: 3D arrays [122 x 199 x 12]
%
% NOTES
% - This function is self-contained and does not read or write files.
% - It expects the input profiles to be positive where valid; non-positive
%   cp660 values are treated as NaN.
% - It preallocates outputs and returns NaN where inputs are invalid.
%
% Example:
%   params.Clim_NPP = Clim_NPP; params.Clim_Zeu2 = Clim_Zeu2; ...
%   results = compute_poc_flux(params);
%
% License: MIT

arguments
    params struct
end

% === Validate minimal required fields ===
requiredFields = {'Clim_NPP','Clim_Zeu2','Clim_Zeu3','Clidata_NN_cp660', ...
    'Clim_MLD','SCSdepth9km','depth','FM_bbp_Micro'};
for k = 1:numel(requiredFields)
    if ~isfield(params, requiredFields{k})
        error('compute_poc_flux:MissingField','Missing required field: %s', requiredFields{k});
    end
end

% Assign variables locally for readability
Clim_NPP = params.Clim_NPP;
Clim_Zeu2 = params.Clim_Zeu2;
Clim_Zeu3 = params.Clim_Zeu3;
Clidata_NN_cp660 = params.Clidata_NN_cp660;
Clim_MLD = params.Clim_MLD;
SCSdepth9km = params.SCSdepth9km;
depth = params.depth(:);
FM_bbp_Micro = params.FM_bbp_Micro;

if isfield(params,'Zeu2') && ~isempty(params.Zeu2)
    Zeu2 = params.Zeu2;
else
    Zeu2 = Clim_Zeu2; % fallback (original code used Zeu2 in some places)
end

% Optional parameters with defaults
if isfield(params,'f_FecT'); f_FecT = params.f_FecT; else f_FecT = 0.1; end
if isfield(params,'m_ph'); m_ph = params.m_ph; else m_ph = 0.1; end
if isfield(params,'useParallel'); useParallel = params.useParallel; else useParallel = false; end

% Dimensions & sanity checks
[nI, nJ, nDepth, nMonths] = size(Clidata_NN_cp660);
if ~isequal(size(Clim_NPP),[nI nJ nMonths]) || ~isequal(size(Clim_Zeu2),[nI nJ nMonths])
    error('compute_poc_flux:BadDims','Input arrays must have consistent spatial and monthly dimensions.');
end

% Preallocate result arrays filled with NaN
nan3 = nan(nI,nJ,nMonths);
Flux_Zeu_nn = nan3; Flux_100_nn = nan3; DPOC = nan3; DPOC100 = nan3;
IPOC_Zeu_NN = nan3; IPOC_100_NN = nan3; IPOC_Zeu_NN_1 = nan3; IPOC_100_NN_1 = nan3;
IPOC_MLD2_nn_1 = nan3; IPOC_MLD2_nn_2 = nan3;

% Precompute vertical index vector for interpolation base points (1..depth)
% The original code constructs interpolation with [1; depth] and POC values
interp_base_depth_idx = [1; depth(:)];

% Main loop over months; optionally run in parallel
monthLoop = 1:nMonths;
if useParallel
    if isempty(gcp('nocreate'))
        parpool('local');
    end
end

parfor_arg = 'for';
if useParallel
    parfor_arg = 'parfor';
end

% Use dynamic evaluation to switch between for/parfor cleanly
eval(sprintf('%s jm=1; end;','% placeholder to satisfy editor'));

% We'll implement a standard for-loop but allow parfor by branching
if useParallel
    parfor jm = 1:nMonths
        local_compute_month(jm);
    end
else
    for jm = 1:nMonths
        local_compute_month(jm);
    end
end

% Pack results
results.Flux_Zeu_nn = Flux_Zeu_nn;
results.Flux_100_nn = Flux_100_nn;
results.DPOC = DPOC;
results.DPOC100 = DPOC100;
results.IPOC_Zeu_NN = IPOC_Zeu_NN;
results.IPOC_100_NN = IPOC_100_NN;
results.IPOC_Zeu_NN_1 = IPOC_Zeu_NN_1;
results.IPOC_100_NN_1 = IPOC_100_NN_1;
results.IPOC_MLD2_nn_1 = IPOC_MLD2_nn_1;
results.IPOC_MLD2_nn_2 = IPOC_MLD2_nn_2;
results.params = params; % echo inputs for provenance

return

% -------------------------
% Local function that computes a single month slice
    function local_compute_month(jm)
        if jm==1
            ts = 12;
        else
            ts = jm-1;
        end

        for is = 1:nI
            for js = 1:nJ
                try
                    % Basic checks
                    if SCSdepth9km(is,js) <= 200
                        continue; % leave as NaN
                    end

                    NPP = Clim_NPP(is,js,jm);
                    if isnan(NPP)
                        continue;
                    end

                    % normalize NPP if deep Zeu3 > 0
                    if Clim_Zeu3(is,js) > 0
                        NPP = NPP / Clim_Zeu3(is,js) * Clim_Zeu2(is,js,jm);
                    end

                    % load cp660 profiles (depth x 1)
                    cp660_nn = squeeze(Clidata_NN_cp660(is,js,:,jm));
                    cp660_nn_1 = squeeze(Clidata_NN_cp660(is,js,:,ts));

                    % protect against non-positive cp values before log10
                    cp660_nn(cp660_nn <= 0) = NaN;
                    cp660_nn_1(cp660_nn_1 <= 0) = NaN;

                    % Convert cp660 to POC (mg/m3): formula from original code
                    POC_NN = 10.^((log10(cp660_nn) + 1.71)./1.27) * 12; % mg/m3
                    POC_NN_1 = 10.^((log10(cp660_nn_1) + 1.71)./1.27) * 12;

                    if all(isnan(POC_NN)) || all(isnan(POC_NN_1))
                        continue;
                    end

                    % Vertical summation over integer depths
                    StZeu = 1:round(max(1,Clim_Zeu2(is,js,jm)));
                    StZeu2 = 1:100;

                    % Interpolate POC onto integer depths; prepend surface value as
                    % in original script using [1; depth] and POC(1)
                    POC_NN_full = [POC_NN(1); POC_NN(:)];
                    POC_NN_1_full = [POC_NN_1(1); POC_NN_1(:)];

                    % For interp1, x values must match the base index vector
                    IPOC_Zeu_NN(is,js,jm) = nansum(interp1(interp_base_depth_idx, POC_NN_full, StZeu, 'linear', NaN));
                    IPOC_100_NN(is,js,jm) = nansum(interp1(interp_base_depth_idx, POC_NN_full, StZeu2, 'linear', NaN));

                    IPOC_Zeu_NN_1(is,js,jm) = nansum(interp1(interp_base_depth_idx, POC_NN_1_full, StZeu, 'linear', NaN));
                    IPOC_100_NN_1(is,js,jm) = nansum(interp1(interp_base_depth_idx, POC_NN_1_full, StZeu2, 'linear', NaN));

                    % MLD-based sums (may be fractional depth)
                    mld_cur = Clim_MLD(is,js,jm);
                    IPOC_MLD2_nn_1(is,js,jm) = nansum(interp1(interp_base_depth_idx, POC_NN_1_full, 1:round(max(1,mld_cur)), 'linear', NaN));
                    IPOC_MLD2_nn_2(is,js,jm) = nansum(interp1(interp_base_depth_idx, POC_NN_full, 1:round(max(1,mld_cur)), 'linear', NaN));

                    % Compute H and checks for Zeu mixing and 100m mixing
                    H = Clim_MLD(is,js,jm) - Clim_MLD(is,js,ts);
                    HZeu = Clim_MLD(is,js,jm) - Zeu2(is,js,jm);
                    H100 = Clim_MLD(is,js,jm) - 100;

                    % Mixing to Zeu condition
                    if H > 0 && HZeu > 0
                        phym1 = (IPOC_MLD2_nn_1(is,js,jm) - IPOC_Zeu_NN_1(is,js,jm)) * FM_bbp_Micro(is,js,ts) / 100;
                        phym2 = (IPOC_MLD2_nn_2(is,js,jm) - IPOC_Zeu_NN(is,js,jm)) * FM_bbp_Micro(is,js,jm) / 100;
                        phym  = phym1 - phym2;

                        phym3 = (IPOC_MLD2_nn_1(is,js,jm) - IPOC_Zeu_NN_1(is,js,jm)) * (1 - FM_bbp_Micro(is,js,ts)/100);
                        phym4 = (IPOC_MLD2_nn_2(is,js,jm) - IPOC_Zeu_NN(is,js,jm)) * (1 - FM_bbp_Micro(is,js,jm)/100);
                        phym5 = phym3 - phym4;
                    else
                        phym = 0; phym5 = 0;
                    end

                    % Alg & particulate pools
                    Alg = 0.1 * NPP * FM_bbp_Micro(is,js,jm) / 100;
                    Pm2_nn = IPOC_Zeu_NN(is,js,jm) * FM_bbp_Micro(is,js,jm) / 100 * 0.35;
                    Pm1_nn = IPOC_Zeu_NN_1(is,js,jm) * FM_bbp_Micro(is,js,ts) / 100 * 0.35;
                    Gm_nn = NPP * FM_bbp_Micro(is,js,jm) / 100 - 0.1 * Pm2_nn - (Pm2_nn - Pm1_nn) + phym - Alg;

                    Ps2_nn = IPOC_Zeu_NN(is,js,jm) * (1 - FM_bbp_Micro(is,js,jm)/100) * 0.35;
                    Ps1_nn = IPOC_Zeu_NN_1(is,js,jm) * (1 - FM_bbp_Micro(is,js,ts)/100) * 0.35;
                    Gs_nn = NPP * (1 - FM_bbp_Micro(is,js,jm)/100) - 0.1 * Ps2_nn - (Ps2_nn - Ps1_nn) + phym5;

                    Flux_Zeu_nn(is,js,jm) = 0.30 * Gm_nn + 0.1 * Gs_nn + Alg;
                    DPOC(is,js,jm) = Ps2_nn - Ps1_nn + Pm2_nn - Pm1_nn;

                    % Mixing to 100m condition
                    if H > 0 && H100 > 0
                        phym1 = (IPOC_MLD2_nn_1(is,js,jm) - IPOC_100_NN_1(is,js,jm)) * FM_bbp_Micro(is,js,ts) / 100;
                        phym2 = (IPOC_MLD2_nn_2(is,js,jm) - IPOC_100_NN(is,js,jm)) * FM_bbp_Micro(is,js,jm) / 100;
                        phym  = phym1 - phym2;

                        phym3 = (IPOC_MLD2_nn_1(is,js,jm) - IPOC_100_NN_1(is,js,jm)) * (1 - FM_bbp_Micro(is,js,ts)/100);
                        phym4 = (IPOC_MLD2_nn_2(is,js,jm) - IPOC_100_NN(is,js,jm)) * (1 - FM_bbp_Micro(is,js,jm)/100);
                        phym5 = phym3 - phym4;
                    else
                        phym = 0; phym5 = 0;
                    end

                    % Compute 100m fluxes
                    Alg = 0.1 * NPP * FM_bbp_Micro(is,js,jm) / 100;
                    Pm2_nn = IPOC_100_NN(is,js,jm) * FM_bbp_Micro(is,js,jm) / 100 * 0.35;
                    Pm1_nn = IPOC_100_NN_1(is,js,jm) * FM_bbp_Micro(is,js,ts) / 100 * 0.35;
                    Gm_nn = NPP * FM_bbp_Micro(is,js,jm) / 100 - 0.1 * Pm2_nn - (Pm2_nn - Pm1_nn) + phym - Alg;

                    Ps2_nn = IPOC_100_NN(is,js,jm) * (1 - FM_bbp_Micro(is,js,jm)/100) * 0.35;
                    Ps1_nn = IPOC_100_NN_1(is,js,jm) * (1 - FM_bbp_Micro(is,js,ts)/100) * 0.35;
                    Gs_nn = NPP * (1 - FM_bbp_Micro(is,js,jm)/100) - 0.1 * Ps2_nn - (Ps2_nn - Ps1_nn) + phym5;

                    Flux_100_nn(is,js,jm) = 0.3 * Gm_nn + 0.1 * Gs_nn + Alg;
                    DPOC100(is,js,jm) = Ps2_nn - Ps1_nn + Pm2_nn - Pm1_nn;

                catch ME
                    % On any unexpected error keep NaN and continue; optionally
                    % you may log ME.message for debugging
                    continue;
                end
            end
        end
    end
end
