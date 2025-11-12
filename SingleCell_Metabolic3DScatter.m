%% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
%
% Script Name: SingleCell_Metabolic3DScatter.m
%
% Description:
%   This script performs single-cell metabolic quantification and 3D
%   visualization based on FTIR single-frequency imaging data.
%   It computes metabolic readouts for each segmented cell and visualizes
%   them in a 3D feature space defined by:
%       x-axis: Lipid uptake      = azide PA / CH2
%       y-axis: Protein synthesis = 13C / (13C + 12C)  (after 2×2 unmixing)
%       z-axis: CH2 / (13C + 12C)
%
%
% Input images (demo filenames shown below):
%   - data_segment.tiff : labeled segmentation mask
%   - data_1616.tif     : 13C amide I channel
%   - data_1650.tif     : 12C amide I channel
%   - data_2094.tif     : azide PA (C–D) channel
%   - data_2914.tif     : CH2 channel
%   (and corresponding files for data2)
%
%  Note:
%     The example datasets (data1 / data2) are demo data for illustrating
%     the workflow. The code itself is production-ready and applicable to
%     experimental single-cell FTIR datasets.
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
%   - MATLAB R2023b or later
%   - Custom function: drawn3D_EE.m  (for 3D error ellipsoid visualization)
%
% Citation:
% If you use this code, please cite the above article.
%
% =========================================================================
%% 
clear all
close all
clc
%% ------------------------------------------------------------------------
% Dataset 1 (data1)
% -------------------------------------------------------------------------
segment1 = imread('data_segment.tiff');
C13_1    = imread('data_1616.tif');
C12_1    = imread('data_1650.tif');
CD_1     = imread('data_2094.tif');
CH2_1    = imread('data_2914.tif');

N1 = max(segment1(:));

for i = 1:N1
    mask = segment1;
    mask(mask ~= i) = 0;
    mask(mask >  0) = 1;

    C13AA_1(i) = sum(sum(C13_1 .* mask)) / 5000;
    C12AA_1(i) = sum(sum(C12_1 .* mask)) / 5000;
    CDsum_1(i) = sum(sum(CD_1  .* mask)) / 5000;
    CH2sum_1(i)= sum(sum(CH2_1 .* mask)) / 5000;
    area1(i)   = sum(sum(mask));
end

C13AA_avg_1 = C13AA_1 ./ area1;
C12AA_avg_1 = C12AA_1 ./ area1;
CDavg_1     = CDsum_1 ./ area1;
CH2avg_1    = CH2sum_1 ./ area1;

% 2×2 spectral unmixing
C12AA_1_unmix = C12AA_avg_1 * 1.281   + C13AA_avg_1 * (-0.548);
C13AA_1_unmix = C12AA_avg_1 * (-0.657)+ C13AA_avg_1 * 1.281;
sum1          = C13AA_1_unmix + C12AA_1_unmix;

ratio_CAA1      = C13AA_1_unmix ./ sum1;
ratio_azidePA1  = CDavg_1      ./ CH2avg_1;
ratio_CH2_1     = CH2avg_1     ./ sum1;

%% ------------------------------------------------------------------------
% Dataset 2 (data2)
% -------------------------------------------------------------------------
segment2 = imread('U251_IDH1_segment.tiff');
C13_2    = imread('U251_IDH1_1616.tif');
C12_2    = imread('U251_IDH1_1650.tif');
CD_2     = imread('U251_IDH1_2094.tif');
CH2_2    = imread('U251_IDH1_2914.tif');

N2 = max(segment2(:));

for i = 1:N2
    mask = segment2;
    mask(mask ~= i) = 0;
    mask(mask >  0) = 1;

    C13AA_2(i) = sum(sum(C13_2 .* mask)) / 5000;
    C12AA_2(i) = sum(sum(C12_2 .* mask)) / 5000;
    CDsum_2(i) = sum(sum(CD_2  .* mask)) / 5000;
    CH2sum_2(i)= sum(sum(CH2_2 .* mask)) / 5000;
    area2(i)   = sum(sum(mask));
end

C13AA_avg_2 = C13AA_2 ./ area2;
C12AA_avg_2 = C12AA_2 ./ area2;
CDavg_2     = CDsum_2 ./ area2;
CH2avg_2    = CH2sum_2 ./ area2;

C12AA_2_unmix = C12AA_avg_2 * 1.281   + C13AA_avg_2 * (-0.548);
C13AA_2_unmix = C12AA_avg_2 * (-0.657)+ C13AA_avg_2 * 1.281;
sum2          = C13AA_2_unmix + C12AA_2_unmix;

ratio_CAA2      = C13AA_2_unmix ./ sum2;
ratio_azidePA2  = CDavg_2      ./ CH2avg_2;
ratio_CH2_2     = CH2avg_2     ./ sum2;

%% ------------------------------------------------------------------------
% Assemble data matrices for ellipsoid calculation
% -------------------------------------------------------------------------
x1 = [ratio_azidePA1', ratio_CAA1', ratio_CH2_1'];
x2 = [ratio_azidePA2', ratio_CAA2', ratio_CH2_2'];

%% ------------------------------------------------------------------------
% Plot 3D scatter + 3D error ellipsoids
% -------------------------------------------------------------------------
colorList = [ ...
    192, 192, 192;   % data1
    103, 168, 205;   % data2
    ] ./ 255;

figure('Position', [100, 100, 700, 600]);

scatter3(ratio_azidePA1, ratio_CAA1, ratio_CH2_1, ...
         25, 'filled', 'MarkerFaceColor', colorList(1,:));
hold on;

scatter3(ratio_azidePA2, ratio_CAA2, ratio_CH2_2, ...
         25, 'filled', 'MarkerFaceColor', colorList(2,:));

% 3D error ellipsoids (requires custom function drawn3D_EE.m)
drawn3D_EE(x1, 'C', 0.9, 'Col', colorList(1,:), 'FA', 0.2, 'PA', false);
drawn3D_EE(x2, 'C', 0.9, 'Col', colorList(2,:), 'FA', 0.2, 'PA', false);

legend({'data1', 'data2'}, ...
       'FontSize', 10, 'FontName', 'Arial', 'Location', 'best');

xlim([0.04 0.12]);
ylim([0.21 0.50]);
zlim([0.05 0.35]);
view(63.5, 24.7);

set(gca, 'FontSize', 14, 'FontName', 'Arial', 'LineWidth', 1);
set(gcf, 'Renderer', 'painters');

xlabel('Lipid uptake (azide PA / CH_2)',                 'FontSize', 14, 'FontName', 'Arial');
ylabel('Protein synthesis (^{13}C / (^{13}C + ^{12}C))', 'FontSize', 14, 'FontName', 'Arial');
zlabel('CH_2 / (^{13}C + ^{12}C)',                       'FontSize', 14, 'FontName', 'Arial');
title('azide PA + ^{13}C DMEM for 24 h',                 'FontSize', 15, 'FontName', 'Arial');

grid on;
box on;



