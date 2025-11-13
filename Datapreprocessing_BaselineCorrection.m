%% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
%
% Script Name: Datapreprocessing_BaselineCorrection.m
% Description:
%   This script demonstrates how to perform rubber-band baseline correction
%   on hyperspectral FTIR imaging data. The workflow includes:
%     1. Loading FTIR data and metadata (Minfo)
%     2. Extracting the desired wavenumber range
%     3. Applying rubber-band baseline correction to each spectrum
%     4. Reshaping the corrected spectra back to the spatial cube
%     5. Updating Minfo and saving the processed result
%
% Input:
%   Example file: 'ExampleData_qualitytest_noisereduction(PCA10).mat'
%     - C     : FTIR hyperspectral cube [Nx, Ny, N_wavenumbers, N_channels]
%     - Minfo : Metadata structure (with fields: File.LWN, File.UWN, File.NOD, File.WVS)
%
% Output:
%   Example output: 'ExampleData_BC.mat'
%     - C     : Baseline-corrected and trimmed FTIR cube
%     - Minfo : Updated metadata corresponding to the trimmed wavenumber region
%
% Notes:
%   - You can modify the spectral trimming indices (d1, d2) according to
%     your desired wavenumber range:
%         1000–1800 cm-1  → d1=1, d2=252
%         2050–2235 cm-1  → d1=1, d2=59
%         2700–3100 cm-1  → d1=1, d2=127
%
% Corresponding paper:
% "VIBRANT-G: Chemical and metabolic profiling of single-cell response to
% genetic perturbations"
%
% Authors: Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang,
%         Xingyuan Bi, Wei Hua, Lixue Shi*
%
% Contact: mqwei23@m.fudan.edu.cn
%
% Dependencies:
%   - MATLAB R2023b (MathWorks, Natick, MA)
%   - Custom function: bc_rubber.m (rubber-band baseline correction)
%
% Citation:
% If you use this code, please cite the above article.
% =========================================================================
%% 
clc;
clear;
close all;

data = load('/path/to/data.mat'); % Example input file
processed = data.C(:,:,:,2); % Select channel 2 (e.g. noise-reduced or PCA component)

Minfo = data.Minfo;
lwn = str2double(data.Minfo.File.LWN);
uwn = str2double(data.Minfo.File.UWN);
nod = str2double(data.Minfo.File.NOD);
wvs = str2double(data.Minfo.File.WVS);
clear data;

[a,b,c] = size(processed);
pspec = reshape(processed, a*b, c);
wave = linspace(lwn, uwn, nod);
%% Trim spectral region and perform rubber-band baseline correction
% Example: 1000–1800 cm^-1 region
d1 = 1;  
d2 = 252;

WN_trim = wave(d1:d2);
spec_trim = pspec(:, d1:d2);
spec_trim_bc = bc_rubber(spec_trim);
[~, e] = size(spec_trim_bc);

C = reshape(spec_trim_bc, a, b, e);
Minfo.File.LWN = num2str(WN_trim(1));
Minfo.File.UWN = num2str(WN_trim(end));
Minfo.File.NOD = num2str(e);

save('/path/to/data_BC.mat', 'C', 'Minfo', '-v7.3');

