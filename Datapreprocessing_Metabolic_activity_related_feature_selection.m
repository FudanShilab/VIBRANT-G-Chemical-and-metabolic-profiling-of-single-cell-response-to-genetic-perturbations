%% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
%
% Script Name: Datapreprocessing_Metabolic_activity_related_feature_selection.m
%
% Description:
%   Performs feature selection on single-cell FTIR spectra by evaluating the
%   Pearson correlation of each wavenumber with two metabolic proxy peaks:
%   13C amide I (1616 cm^-1) and azido (2094 cm^-1). A composite score is
%   derived from these correlations and used to rank spectral features.
%   The top 100 ranked wavenumbers are exported for downstream analysis.
%
% Data format (processed_data.csv):
%   - Row 1: wavenumber axis (cm^-1), from column 2 onward
%   - Col 1: cell ID (integer index or label)
%   - Col 2..end: single-cell spectra (rows = cells, columns = wavenumbers)
%
% Corresponding paper:
%   "VIBRANT-G: Chemical and metabolic profiling of single-cell response
%    to genetic perturbations"
%   Authors: Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang,
%            Xingyuan Bi, Wei Hua, Lixue Shi*
%
% Contact: mqwei23@m.fudan.edu.cn
%
% Dependencies:
%   - MATLAB R2023b (MathWorks, Natick, MA)
%   - Built-in corr() function for Pearson correlation
%
% Citation:
%   If you use this code, please cite the above article.
%
% =========================================================================
%% 
clear; close all; clc;
%% ========================== USER SETTINGS ================================
in_csv       = '/path/to/processed_data.csv';
anchor1_idx  = 195;    % Index of 13C amide I
anchor2_idx  = 265;    % Index of azido PA
TopN         = 100;    % Number of top-ranked features to export
TopK_report  = 100;    % Number of top features to include in ranking table
out_dir      = pwd;    % Output directory (default: current folder)

%% ========================== LOAD & PARSE DATA ============================
if ~isfile(in_csv)
    error('Input file not found: %s', in_csv);
end

raw = readmatrix(in_csv);
if isempty(raw)
    error('Input CSV is empty or unreadable.');
end

wavenumbers = raw(1, 3:end);      % 1 x F
X           = raw(2:end, 3:end);  % N x F (cells × features)
[N, F] = size(X);
fprintf('Loaded %d cells × %d features.\n', N, F);

%% =================== REMOVE INVALID OR ZERO-VARIANCE FEATURES ============
valid_cols = all(isfinite(X), 1);
X = X(:, valid_cols);
wavenumbers = wavenumbers(valid_cols);
F = size(X, 2);

std_feature = std(X, 0, 1);
nonzero_mask = std_feature > 1e-12;
X = X(:, nonzero_mask);
wavenumbers = wavenumbers(nonzero_mask);
F = size(X, 2);
fprintf('After filtering: %d valid features remain.\n', F);

if any(~ismember([anchor1_idx, anchor2_idx], 1:F))
    error('Anchor indices are out of range (1–%d after filtering).', F);
end

%% ==================== FEATURE–FEATURE CORRELATION ========================
C = corr(X, 'type', 'Pearson', 'rows', 'pairwise');  % F × F correlation matrix
a1 = abs(C(:, anchor1_idx));                         % Correlation with 13C amide I
a2 = abs(C(:, anchor2_idx));                         % Correlation with azido PA

%% =========================== SCORING & RANKING ==========================
score = a1 + a2 - (a1 .* a2);                        % Composite score
[score_sorted, sort_idx] = sort(score, 'descend');    % Sort descending

%% ============================ EXPORT RESULTS =============================
TopN = min(TopN, F);
sel_by_score = sort_idx(1:TopN);
sel_in_order = sort(sel_by_score);

Wsel = wavenumbers(sel_in_order);
Xsel = X(:, sel_in_order);
out_matrix = [Wsel; Xsel];

writematrix(out_matrix, fullfile(out_dir, 'top_features_matrix.csv'));

K = min(TopK_report, F);
T = table( ...
    sort_idx(1:K), ...
    wavenumbers(sort_idx(1:K)).', ...
    score_sorted(1:K), ...
    'VariableNames', {'feature_index','wavenumber_cm_1','score'} );

writetable(T, fullfile(out_dir, 'top_feature_indices.csv'));

fprintf('Export complete:\n');
fprintf('  - top_features_matrix.csv\n');
fprintf('  - top_feature_indices.csv\n');

%% ========================== VISUALIZATION ================================
avg = mean(X, 1);
if exist('spcnormalize', 'file')
    avg = spcnormalize(avg);
else
    % Fallback: simple min–max normalization
    avg = (avg - min(avg)) / (max(avg) - min(avg) + eps);
end

figure('Color','w','Position',[100 500 1100 300]);
scatter(1:F, avg, 40, score, 'filled');
cb = colorbar;
colormap(turbo);  % MATLAB built-in high-contrast colormap
xlabel('Feature Index (1 to F)', 'FontSize', 11);
ylabel('Mean Intensity', 'FontSize', 11);
title('Feature Score Scatter Plot (Color = Feature Score)', 'FontSize', 12);
ylabel(cb, 'Feature Score');
grid on; box on;

fprintf('Visualization complete.\n');
