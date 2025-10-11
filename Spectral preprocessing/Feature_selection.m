%% =========================================================================
% Script Name: Metabolic activity related feature selection.m
% Description:
%   This script performs feature selection on single-cell FTIR spectra by
%   evaluating the Pearson correlation of each wavenumber with two key
%   metabolic proxies (13C amide I at 1616 cm^-1 and azido at 2096 cm^-1).
%   A composite score is derived from these correlations, and features are
%   ranked accordingly. The top 100 ranked wavenumbers are exported for
%   downstream analysis.

% Data format (processed_data.csv):
%   - Row 1: wavenumber axis (cm^-1), from column 2 onward
%   - Col 1: cell ID (integer index or label)
%   - Col 2..end: single-cell spectra (rows = cells, columns = wavenumbers)

% Corresponding paper:
%   "VIBRANT-G: Chemical and metabolic profiling of single-cell response
%    to genetic perturbations"
%   Authors: Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang,
%            Dan Ye, Wei Hua, Lixue Shi*
%
% Dependencies:
%   - MATLAB R2023b (MathWorks, Natick, MA)
%   - Pearson correlation computed using MATLAB built-in corr() function
%
% Citation:
%   If you use this code, please cite the above article.
%
% License:For academic use only.
% =========================================================================
%% 
clear; close all; clc;
%% ===== USER SETTINGS (edit here) =====
in_csv       = "/path/to/your/data.csv"; 
anchor1_idx  = 195;    % Index of 13C amide I 
anchor2_idx  = 265;    % Index of azido PA
TopN         = 100;    % Number of top-ranked features to export 
TopK_report  = 100;    % Number of top features to report in the ranking table (by score order)

%% ===== LOAD & BASIC PARSE =====
raw = readmatrix(in_csv);
wavenumbers = raw(1, 3:end);           % 1 x F
X           = raw(2:end, 3:end);       % N x F (cells x features)

[N, F] = size(X);
if any(~isfinite(wavenumbers), 'all') || any(~isfinite(X), 'all')
    error('Non-finite values in input CSV.');
end
if any(~ismember([anchor1_idx, anchor2_idx], 1:F))
    error('Anchor indices out of range (must be within 1..F).');
end

%% ===== FEATURE–FEATURE CORRELATION =====
C = corr(X);            % F x F
a1 = abs(C(:, anchor1_idx));   
a2 = abs(C(:, anchor2_idx));  

%% ===== SCORE & RANK =====
score = a1 + a2 - (a1 .* a2);        
[score_sorted, sort_idx] = sort(score, 'descend');  

%% ===== SELECT & EXPORT =====
% (1) TopN kept in original spectral order (ascending by column index)
TopN = min(TopN, F);
sel_by_score = sort_idx(1:TopN);          
sel_in_order = sort(sel_by_score);        

Wsel = wavenumbers(sel_in_order);      
Xsel = X(:, sel_in_order);           
out_matrix = [Wsel; Xsel];             
writematrix(out_matrix, 'top_features_matrix.csv');

% (2) Ranking table for reporting (TopK_report by score order)
K = min(TopK_report, F);
T = table( ...
    sort_idx(1:K), ...
    wavenumbers(sort_idx(1:K)).', ...
    score_sorted(1:K), ...
    'VariableNames', {'feature_index','wavenumber_cm_1','score'} );
writetable(T, 'top_feature_indices.csv');

fprintf('Done. Wrote:\n  - top_features_matrix.csv\n  - top_feature_indices.csv\n');
