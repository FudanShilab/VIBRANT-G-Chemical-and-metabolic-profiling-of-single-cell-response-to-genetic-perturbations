%% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
%
% Script Name:spcnormalize.m
% Description: Performs normalization
%
% Corresponding paper:
% "VIBRANT-G: Chemical and metabolic profiling of single-cell response to
%             genetic perturbations"

% Authors:  Minqian Wei, Minjie Fu, Tongqi Wang, Yuchen Sun, Mingyu Wang
%           Xingyuan Bi, Wei Hua, Lixue Shi*

% Contact: mqwei23@m.fudan.edu.cn
%
% Dependencies:
%   - MATLAB R2023b or later
%
% Citation:
% If you use this code, please cite the above article.
% =========================================================================
%% 
function [ spectra ] = spcnormalize(spectrainput)
spectra = spectrainput;
spectra = spectra - min(spectra,[],2);
spectra = spectra ./ max(spectra,[],2);
end



