function [h, info] = drawn3D_EE(data, varargin)
% =========================================================================
% Copyright (c) 2025 Lixue Shi and collaborators
% Licensed under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0)
%
% Script Name:drawn3D_EE
%
% Description:
%   Convenience wrapper for drawing 3D error ellipsoids from single-cell
%   feature data. This function provides a compact interface with short
%   option names and forwards the corresponding full-parameter set to
%   an internal function drawErrorEllipsoid3D().
%
%   Typical use (3D data, N × 3):
%       drawn3D_EE(X, 'C', 0.99, 'Col', [0.9 0.2 0.2], 'FA', 0.25);
%
%   where:
%       X   : N × 3 matrix, each row is one data point in 3D space.
%
%   Shorthand options (mapped internally to full parameter names):
%       'C'      -> 'Confidence'        (e.g., 0.9, 0.95, 0.99)
%       'Col'    -> 'Color'             (RGB triplet, e.g., [0.9 0.2 0.2])
%       'FA'     -> 'FaceAlpha'         (0–1, ellipsoid face transparency)
%       'EC'     -> 'EdgeColor'
%       'EA'     -> 'EdgeAlpha'
%       'LW'     -> 'LineWidth'
%       'Res'    -> 'Resolution'        (mesh density)
%       'Ax'     -> 'Axes'
%       'Eq'     -> 'EqualAxes'
%       'Axes'   -> 'Axes'
%       'PA'     -> 'DrawPrincipalAxes' (true/false)
%
%   This wrapper always sets 'EqualAxes' to false by default so that the
%   ellipsoid is rendered in the native scaling of the plotted features.
%
% Inputs:
%   data      : N × 3 numeric array, rows are observations, columns are
%               features (e.g., [ratio_azidePA, ratio_CAA, ratio_CH2]).
%   varargin  : name–value pairs using either shorthand (e.g. 'C','Col',
%               'FA', 'EC', 'EA', 'LW', 'Res', 'Ax', 'Eq', 'Axes', 'PA')
%               or the full parameter names accepted by drawErrorEllipsoid3D.
%
% Outputs:
%   h    : Graphics handle(s) returned by drawErrorEllipsoid3D.
%   info : Struct with ellipsoid parameters (center, radii, rotation, cov).
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
%   - MATLAB R2023b or later
%
% Citation:
% If you use this code, please cite the above article.
% =========================================================================

    % Shorthand → full-name mapping
    map = containers.Map( ...
        {'C','Col','FA','EC','EA','LW','Res','Ax','Eq','Axes','PA'}, ...
        {'Confidence','Color','FaceAlpha','EdgeColor','EdgeAlpha', ...
         'LineWidth','Resolution','Axes','EqualAxes','Axes', ...
         'DrawPrincipalAxes'} ...
    );

    % Expand shorthand options into full parameter names
    args = {};
    k = 1;
    while k <= numel(varargin)
        if ischar(varargin{k}) || isstring(varargin{k})
            key = char(varargin{k});
            if isKey(map, key)
                % Replace shorthand with full name
                args(end+1:end+2) = {map(key), varargin{k+1}}; %#ok<AGROW>
                k = k + 2;
                continue;
            end
        end
        % Pass through unmodified if not a recognized shorthand
        args{end+1} = varargin{k}; %#ok<AGROW>
        k = k + 1;
    end

    % Delegate to internal ellipsoid-drawing function.
    % Force EqualAxes = false by default (unless overridden in args).
    [h, info] = drawErrorEllipsoid3D(data, args{:}, 'EqualAxes', false);
end


%% ========================================================================
% Internal function: drawErrorEllipsoid3D
% ========================================================================
function [h, info] = drawErrorEllipsoid3D(data, varargin)
% drawErrorEllipsoid3D  Draw a confidence ellipsoid (error ellipsoid) in 3D.
%
% USAGE:
%   [h, info] = drawErrorEllipsoid3D(data, 'Confidence', 0.95, ...
%                                    'Color', [0.2 0.6 0.9], ...);
%
% Required input:
%   data : N × 3 matrix (each row one point). Rows containing NaN are
%          automatically removed.
%
% Name–Value optional parameters:
%   'Confidence'       : Confidence level in (0,1). Default: 0.95
%                        Internally uses chi2inv(Confidence, 3).
%   'Color'            : Ellipsoid face color (RGB), default [0.2 0.6 0.9].
%   'FaceAlpha'        : Face transparency in [0,1], default 0.2.
%   'EdgeColor'        : Edge color, default 'none'.
%   'EdgeAlpha'        : Edge transparency in [0,1], default 0.8.
%   'LineWidth'        : Edge line width, default 1.
%   'Resolution'       : Sphere mesh resolution (integer ≥ 6), default 50.
%   'Axes'             : Axes handle to plot into, default gca.
%   'EqualAxes'        : Logical; if true, calls axis equal vis3d. Default true.
%   'DrawPrincipalAxes': Logical; if true, draw principal axes lines. Default false.
%
% Outputs:
%   h    : Handle to the surf object representing the ellipsoid.
%   info : Struct with fields:
%            - center : 1 × 3 mean vector
%            - radii  : 1 × 3 radii along principal axes
%            - R      : 3 × 3 rotation matrix (columns = principal directions)
%            - C      : 3 × 3 covariance matrix
%
% Note:
%   This function is intended as the core ellipsoid renderer. For a more
%   compact user-facing interface with shorthand options, use drawn3D_EE().

    % --- Parse inputs ---
    p = inputParser;
    p.addParameter('Confidence', 0.95, @(x) isscalar(x) && x > 0 && x < 1);
    p.addParameter('Color', [0.2 0.6 0.9], @(x) isnumeric(x) && numel(x) == 3);
    p.addParameter('FaceAlpha', 0.2, @(x) isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('EdgeColor', 'none');
    p.addParameter('EdgeAlpha', 0.8, @(x) isscalar(x) && x >= 0 && x <= 1);
    p.addParameter('LineWidth', 1, @(x) isscalar(x) && x >= 0);
    p.addParameter('Resolution', 50, @(x) isscalar(x) && x >= 6);
    p.addParameter('Axes', gca, @(x) ishandle(x) && strcmp(get(x, 'Type'), 'axes'));
    p.addParameter('EqualAxes', true, @(x) islogical(x));
    p.addParameter('DrawPrincipalAxes', false, @(x) islogical(x));
    p.parse(varargin{:});
    S = p.Results;

    % --- Data checks and cleaning ---
    if size(data, 2) ~= 3
        error('drawErrorEllipsoid3D:InvalidData', 'Input "data" must be an N×3 matrix.');
    end
    data = data(~any(isnan(data), 2), :);
    if size(data, 1) < 2
        error('drawErrorEllipsoid3D:NotEnoughData', 'Not enough valid data points (<2).');
    end

    % --- Mean & covariance ---
    center = mean(data, 1, 'omitnan');
    C = cov(data, 'omitrows');
    if rank(C) < 3
        % Light regularization to avoid singular covariance matrices
        C = C + 1e-12 * eye(3);
    end

    % --- Eigen decomposition (ensure symmetry) ---
    [R, D] = eig((C + C') / 2);
    evals = max(diag(D), 0);  % Clamp negative numerical noise to zero

    % --- Chi-square scaling (d = 3) ---
    chiVal = chi2inv(S.Confidence, 3);
    radii  = sqrt(chiVal .* evals(:)');  % 1×3

    % --- Sphere → ellipsoid mapping ---
    [xs, ys, zz] = sphere(S.Resolution);
    Sx = R * diag(radii);
    P  = [xs(:), ys(:), zz(:)] * Sx.';  % Transform unit sphere
    Xs = reshape(P(:,1), size(xs)) + center(1);
    Ys = reshape(P(:,2), size(ys)) + center(2);
    Zs = reshape(P(:,3), size(zz)) + center(3);

    % --- Plotting ---
    ax = S.Axes;
    axes(ax); %#ok<LAXES>
    h = surf(ax, Xs, Ys, Zs, ...
        'FaceColor',  S.Color, ...
        'FaceAlpha',  S.FaceAlpha, ...
        'EdgeColor',  S.EdgeColor, ...
        'EdgeAlpha',  S.EdgeAlpha, ...
        'LineWidth',  S.LineWidth);
    hold(ax, 'on');

    % Optional: draw principal axes
    if S.DrawPrincipalAxes
        for k = 1:3
            v = R(:, k) * radii(k);
            plot3(ax, ...
                  [center(1) - v(1), center(1) + v(1)], ...
                  [center(2) - v(2), center(2) + v(2)], ...
                  [center(3) - v(3), center(3) + v(3)], ...
                  '-', 'Color', S.Color, 'LineWidth', 1.5, ...
                  'HandleVisibility', 'off');
        end
    end

    if S.EqualAxes
        axis(ax, 'equal');
        axis(ax, 'vis3d');
    end

    % --- Output info struct ---
    info.center = center;
    info.radii  = radii;   % [a b c]
    info.R      = R;       % principal directions
    info.C      = C;       % covariance matrix
end
