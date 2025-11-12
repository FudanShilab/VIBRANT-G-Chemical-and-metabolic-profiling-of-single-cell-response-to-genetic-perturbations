
%% =========================================================================
% Script Name: spectral shift correction.m

% Description:
%   Penalized Reference Matching (PRM) for segmented spectral shift 
%   correction on single-cell FTIR data. All inputs are CSVs with:
%     - row = single cell (one spectrum per row)
%     - col1 = cell_id (string)
%     - col2..end = intensity at wavenumbers specified by header row
%   The script:
%     1) Reads reference batch (e.g., WT) and target batch (to be aligned)
%     2) Unifies to a 0.5 cm^-1 grid via spline interpolation
%     3) Computes mean spectra of ref/target
%     4) For each segment (FP / TB / CH), searches candidate shifts Δx
%        to maximize:  S = cos(u, v) - alpha * (Δx)^2
%     5) Applies the 3 optimal segmental shifts to the entire target batch
%     6) Exports: optimal shifts, aligned mean spectrum, per-cell aligned
%        spectra
%
% Corresponding paper:
% "VIBRANT-G: Chemical and metabolic profiling of single-cell response to
% genetic perturbations"
%
% Authors: Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang,
%          Dan Ye, Wei Hua, Lixue Shi*

% Contact: mqwei23@m.fudan.edu.cn
%
% Dependencies: MATLAB R2023b+. 
%               Using the built-in interp1 function in MATLAB
% License: For academic use only.
% =========================================================================
%% 
clear; clc; close all;
%% =========================
%           CONFIG
% =========================
ref_csv = 'ref_spectrum.csv';        % reference batch (e.g., WT)
tgt_csv = 'target_spectrum.csv';     % target batch to align

dw = 0.5;                         % resampling step (cm^-1)
interp_kind = 'spline';           
extrap_val  = 'extrap';           

seg_ranges = [1000 1800;          % Fingerprint Region(FP)
              2050 2235;          % Triple-bond Region(TB)
              2700 3100];         % CH stretch Region(CH)

dx_min  = -15;                    
dx_max  =  15;                    
dx_step = 0.5;                    

alpha = 1e-4;                     
do_l2_normalize = true;           

outdir = 'prm_outputs';
save_aligned_batch   = true;      
save_aligned_mean    = true;      
save_optimal_shifts  = true;      

%% =========================
%        I/O & PARSING
% =========================
fprintf('Reading CSVs...\n');
refT = readtable(ref_csv, 'TextType', 'string');
tgtT = readtable(tgt_csv, 'TextType', 'string');

w_ref = str2double(refT.Properties.VariableNames(2:end));
w_tgt = str2double(tgtT.Properties.VariableNames(2:end));

Xref = table2array(refT(:,2:end));   
Xtgt = table2array(tgtT(:,2:end));   

cell_id_tgt = tgtT{:,1};

ref_raw = mean(Xref, 1).';
tgt_raw = mean(Xtgt, 1).';

%% =========================
%     UNIFY TO 0.5 CM-1
% =========================
lo = max(min(w_ref), min(w_tgt));
hi = min(max(w_ref), max(w_tgt));
wg = (ceil(lo/dw)*dw) : dw : (floor(hi/dw)*dw);  

ref_mean = interp1(w_ref, ref_raw, wg, interp_kind, extrap_val);
tgt_mean = interp1(w_tgt, tgt_raw, wg, interp_kind, extrap_val);

if do_l2_normalize
    ref_mean = ref_mean ./ max(norm(ref_mean), 1e-12);
    tgt_mean = tgt_mean ./ max(norm(tgt_mean), 1e-12);
end

Xtgt_res = interp1(w_tgt, Xtgt.', wg, interp_kind, extrap_val).';  
if do_l2_normalize
    Xtgt_res = Xtgt_res ./ max(vecnorm(Xtgt_res, 2, 2), 1e-12);
end

%% =========================
%   PRM: SEARCH OPTIMAL Δx
% =========================
fprintf('Searching optimal shifts per segment...\n');
dx_candidates = dx_min:dx_step:dx_max;
best_dx = zeros(size(seg_ranges,1),1);

for s = 1:size(seg_ranges,1)
    idx = (wg >= seg_ranges(s,1)) & (wg <= seg_ranges(s,2));
    v = ref_mean(idx);
    if do_l2_normalize, v = v ./ max(norm(v),1e-12); end

    bestScore = -Inf; bestShift = 0;
    for dx = dx_candidates
        u = interp1(wg - dx, tgt_mean, wg(idx), interp_kind, extrap_val);
        if do_l2_normalize, u = u ./ max(norm(u),1e-12); end

        score  = dot(u, v) - alpha * (dx^2);
        if score > bestScore
            bestScore = score;
            bestShift = dx;
        end
    end
    best_dx(s) = bestShift;
end

fprintf('Optimal Δx (cm^-1): [%.2f %.2f %.2f]\n', best_dx);

%% =========================
%  APPLY Δx TO ALL TARGETS
% =========================
Xtgt_aligned = Xtgt_res;
for s = 1:size(seg_ranges,1)
    idx = (wg >= seg_ranges(s,1)) & (wg <= seg_ranges(s,2));
    dx  = best_dx(s);
    Xtgt_aligned(:, idx) = interp1(wg - dx, Xtgt_res.', wg(idx), interp_kind, extrap_val).';
end
if do_l2_normalize
    Xtgt_aligned = Xtgt_aligned ./ max(vecnorm(Xtgt_aligned, 2, 2), 1e-12);
end

tgt_aligned_mean = mean(Xtgt_aligned, 1).';
if do_l2_normalize
    tgt_aligned_mean = tgt_aligned_mean ./ max(norm(tgt_aligned_mean), 1e-12);
end

%% =========================
%         OUTPUTS
% =========================
if ~exist(outdir,'dir'), mkdir(outdir); end

if save_optimal_shifts
    Tdx = table(["FP";"TB";"CH"], best_dx, ...
        'VariableNames', {'segment','optimal_shift_cm_1'});
    writetable(Tdx, fullfile(outdir, 'optimal_shifts.csv'));
end

if save_aligned_mean
    Tout = table(wg(:), ref_mean(:), tgt_mean(:), tgt_aligned_mean(:), ...
        'VariableNames', {'wavenumber','ref_mean','tgt_mean','tgt_aligned_mean'});
    writetable(Tout, fullfile(outdir, 'aligned_mean.csv'));
end

if save_aligned_batch
    Tall = array2table(Xtgt_aligned, 'VariableNames', cellstr(string(wg)));
    Tall = addvars(Tall, string(cell_id_tgt), 'Before', 1, 'NewVariableNames', 'cell_id');
    writetable(Tall, fullfile(outdir, 'aligned_batch.csv'));
end

fprintf('\nDone. Outputs saved in: %s\n', outdir);
