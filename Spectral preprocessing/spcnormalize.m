%% =========================================================================
% Script Name:spcnormalize.m
% Description: Performs normalization
%
% Corresponding paper:
% "VIBRANT-G: Chemical and metabolic profiling of single-cell response to
%             genetic perturbations"

% Authors:  Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang
%           Dan Ye, Wei Hua, Lixue Shi*

% Email:mqwei23@m.fudan.edu.cn
%
% Dependencies:
%   - MATLAB R2024a or later
%
% Citation:
% If you use this code, please cite the above article.
%
% License: For academic use only.
% =========================================================================
function [ spectra ] = spcnormalize(spectrainput)
spectra = spectrainput;
spectra = spectra - min(spectra,[],2);
spectra = spectra ./ max(spectra,[],2);
end



