function ax = plot_cross_section(M, opts)
% PLOT_CROSS_SECTION
% -------------------------------------------------------------------------
% Plots a single cross-section of an unstructured hexahedral mesh by
% intersecting it with ONE user-defined plane. The slice preserves the
% original cell geometry by computing polygonal intersections of the
% selected plane with each hexahedral cell.
%
% The function supports plotting resistivity or conductivity stored at
% nodes and can produce a 2D projected view or a 3D view.
%
% IMPORTANT
%   - Exactly ONE plane must be specified per call using:
%       'x', 'y', 'z', or 'point'.
%   - To plot multiple planes, call this function multiple times.
%
% NAME-VALUE OPTIONS
%   'view'     : '2D' or '3D' (default: '2D')
%   'property' : 'conductivity'|'cond'|'sigma' or 'resistivity'|'res'|'rho'
%                (default: 'resistivity')
%   'axes'     : Axes handle. If empty, the function creates an axes.
%
% PLANE KEYS (exactly one must be provided)
%   'x'     : scalar plane value for x = const
%   'y'     : scalar plane value for y = const
%   'z'     : scalar plane value for z = const
%   'point' : [A;B] with size 2x2 or 2x3
%
% ASSUMPTIONS / LIMITATIONS
%   - Only VTK_HEXAHEDRON cells (type 12) are used for slicing.
%   - Intersections are computed with linear interpolation along edges.
%
% SYNTAX
%   plot_cross_section(M, 'x',      0)
%   plot_cross_section(M, 'view',   '2D', 'property', 'rho', 'y', 10)
%   plot_cross_section(M, 'view',   '3D', 'property', 'cond', 'z', -2)
%   plot_cross_section(M, 'property', 'res', 'point', [-10 -10; 10 10])
%   plot_cross_section(M, 'view',   '2D', 'property', 'rho', 'x', 0, 'axes', ax)
% -------------------------------------------------------------------------
% UPDATE LOG
%   version         description
%   1.0             First release
% -------------------------------------------------------------------------
%   Author:         César Daniel Castro
%   Affiliation:    Institute of Geophysics, Czech Academy of Science
%   First release:  08-Dec-2025
%   Last update:    08-Dec-2025
% -------------------------------------------------------------------------

    arguments
        M struct

        % Core options
        opts.view       (1,1) string {mustBeMember(opts.view,       ["2D", "3D", "2d", "3d"])} = "2D"
        opts.property   (1,1) string {mustBeMember(opts.property,   ["conductivity", "cond", "sigma", ...
                                                                     "resistivity", "res", "rho"])} = "resistivity"
        opts.axes       = []

        % Plane keys
        opts.x          (1,1) double = NaN
        opts.y          (1,1) double = NaN
        opts.z          (1,1) double = NaN
        opts.point      double   = NaN
    end

    % Validate axes handle
    % ---------------------------------------------------------------------
    AX = opts.axes;
    if ~isempty(AX) && ~isgraphics(AX, 'axes')
        error('''axes'' must be a valid axes handle or empty [].');
    end

    % Normalize tokens
    % ---------------------------------------------------------------------
    view = upper(string(opts.view));
    prop = lower(string(opts.property));

    % Determine requested plane (exactly one)
    % ---------------------------------------------------------------------
    hasX = ~isnan(opts.x);
    hasY = ~isnan(opts.y);
    hasZ = ~isnan(opts.z);
    hasPoint = ~(isnumeric(opts.point) && isscalar(opts.point) && isnan(opts.point));
    nPlanes = hasX + hasY + hasZ + hasPoint;

    if nPlanes ~= 1
        error('Exactly ONE plane must be specified per call using ''x'', ''y'', ''z'', or ''point''.');
    end

    if hasX
        plane = parse_plane('x', opts.x);
    elseif hasY
        plane = parse_plane('y', opts.y);
    elseif hasZ
        plane = parse_plane('z', opts.z);
    else
        AB = opts.point;
        if ~(isnumeric(AB) && size(AB,1) == 2 && (size(AB,2) == 2 || size(AB,2) == 3))
            error('''point'' must be a 2x2 or 2x3 numeric array [A; B].');
        end
        plane = parse_plane('point', AB);
    end

    % Extract mesh
    % ---------------------------------------------------------------------
    P = M.points;       % [Npts x 3]
    C = M.cells;        % [Ncells x 8]
    types = M.cellTypes;

    if isempty(P) || isempty(C) || isempty(types)
        error('Mesh struct M is missing required fields (.points, .cells, .cellTypes).');
    end

    globalScale = max(range(P, 1));
    if ~isfinite(globalScale) || globalScale <= 0
        globalScale = 1;
    end

    % Select property
    % ---------------------------------------------------------------------
    if ~isfield(M, 'pointData')
        error('Struct M does not contain .pointData.');
    end

    switch prop
        case {"conductivity","cond","sigma"}
            if ~isfield(M.pointData, 'conductivity')
                error('M.pointData does not contain ''conductivity''.');
            end
            values = M.pointData.conductivity;
            propLabel = 'log_{10} Conductivity [S/m]';

        case {"resistivity","res","rho"}
            if ~isfield(M.pointData, 'resistivity')
                error('M.pointData does not contain ''resistivity''.');
            end
            values = M.pointData.resistivity;
            propLabel = 'log_{10} Resistivity [\Omegam]';
    end

    % Filter to hex cells only
    % ---------------------------------------------------------------------
    isHex = (types == 12);
    C = C(isHex, :);

    if isempty(C)
        error('No VTK_HEXAHEDRON (type 12) cells were found in this mesh.');
    end

    % Hex edge list (local indices)
    % ---------------------------------------------------------------------
    hexEdges = [
                    1 2; 2 3; 3 4; 4 1;   % bottom face
                    5 6; 6 7; 7 8; 8 5;   % top face
                    1 5; 2 6; 3 7; 4 8    % vertical edges
                ];

    % Dispatch by view
    % ---------------------------------------------------------------------
    switch view
        case "3D"
            if ~isempty(AX)
                ax = AX;
                hold(ax, 'on');
            else
                fig = figure('Name', plane.label, ...
                             'WindowStyle', 'docked');
                clf(fig);
                ax = axes('Parent', fig);
                hold(ax, 'on');
            end

            cvals_plane_3d = draw_plane_3d(ax, ...
                                           P, C, values, ...
                                           plane.P0, ...
                                           plane.n, ...
                                           hexEdges, ...
                                           plane.proj, ...
                                           plane.value, ...
                                           globalScale);
            format_3d_axes(ax);

            if ~isempty(cvals_plane_3d)
                try
                    clim(ax, [min(cvals_plane_3d) max(cvals_plane_3d)]);
                catch
                end
            end

            % Colorbar handling
            cb = colorbar(ax);
            cb.Label.String = propLabel;

        case "2D"
            if ~isempty(AX)
                ax = AX;
                hold(ax, 'on');
            else
                fig = figure('Name', plane.label, ...
                             'WindowStyle', 'docked');
                clf(fig);
                ax = axes('Parent', fig);
                hold(ax, 'on');
            end

            cvals_plane_2d = draw_plane_2d(ax, ...
                                           P, C, values, ...
                                           plane.P0, ...
                                           plane.n, ...
                                           hexEdges, ...
                                           plane.proj, ...
                                           plane.value, ...
                                           plane, ...
                                           globalScale);

            format_2d_axes(ax, plane.proj, plane.label);
            if strcmp(plane.proj, 'point') && isfield(plane, 'L')
                xlim(ax, [0 plane.L]);
            end

            if ~isempty(cvals_plane_2d)
                try
                    clim(ax, [min(cvals_plane_2d) max(cvals_plane_2d)]);
                catch
                end
            end

            cb = colorbar(ax);
            cb.Label.String = propLabel;
    end

    try
        colormap(ax, flipud(colors('spectral')));
    catch
        colormap(ax, 'parula');
    end

end

function plane = parse_plane(key, val)
% PARSE_PLANE
% -------------------------------------------------------------------------
% Builds a plane definition struct for a single plane key.
% -------------------------------------------------------------------------

    key = lower(string(key));

    switch key
        case "x"
            x0 = val;
            plane.P0 = [x0, 0, 0];
            plane.n  = [1, 0, 0];
            plane.proj = 'x';
            plane.value = x0;
            plane.label = sprintf('Cross-section X = %.2g [km]', x0);

        case "y"
            y0 = val;
            plane.P0 = [0, y0, 0];
            plane.n  = [0, 1, 0];
            plane.proj = 'y';
            plane.value = y0;
            plane.label = sprintf('Cross-section Y = %.2g [km]', y0);

        case "z"
            z0 = val;
            plane.P0 = [0, 0, z0];
            plane.n  = [0, 0, 1];
            plane.proj = 'z';
            plane.value = z0;
            plane.label = sprintf('Plane view Z = %.2g [km]', z0);

        case {'point','profile','line'}
            AB = val;
        
            if size(AB,1) ~= 2
                error('Value for ''point'' must be [A; B] with 2 rows.');
            end
        
            if size(AB,2) == 2
                A = [AB(1,:), 0];
                B = [AB(2,:), 0];
            elseif size(AB,2) == 3
                A = AB(1,:);
                B = AB(2,:);
            else
                error('Value for ''point'' must be 2x2 or 2x3.');
            end
        
            d = B - A;
            if norm(d) < 1e-12
                error('Points A and B for ''point'' must be distinct.');
            end
        
            % Use only horizontal direction to build a vertical plane
            d_xy = [d(1), d(2), 0];
            if norm(d_xy) < 1e-12
                d_xy = [1, 0, 0];
            end
        
            % Plane normal: cross(profile direction, vertical)
            v  = [0, 0, 1];
            n  = cross(d_xy, v);
            if norm(n) < 1e-12
                n = [0, 1, 0];
            end
        
            % Profile axis (unit) and length in map view
            eProfile = d_xy / norm(d_xy);            % unit vector along A->B (horizontal)
            L        = norm(d_xy(1:2));              % segment length (XY)
        
            P0    = [A(1), A(2), 0];
            label = sprintf('Profile: A(%.3g,%.3g) - B(%.3g,%.3g)', ...
                            A(1), A(2), B(1), B(2));
        
            pval = NaN;

            plane.P0    = P0;
            plane.n     = n(:).' / norm(n);
            plane.proj  = 'point';
            plane.value = pval;
            plane.label = label;
            
            % Extra geometry needed for correct distance axis
            plane.A        = A;
            plane.B        = B;
            plane.eProfile = eProfile;
            plane.L        = L;

        otherwise
            error('Unknown plane key "%s".', key);
    end

    % Normalize normal
    plane.n = plane.n(:).';
    plane.n = plane.n / norm(plane.n);

end

function format_3d_axes(ax)
% FORMAT_3D_AXES
% -------------------------------------------------------------------------
    set(ax, 'FontName', 'Verdana', 'FontSize', 12);
    axis(ax, 'equal', 'tight');
    xlabel(ax, 'Easting [km]');
    ylabel(ax, 'Northing [km]');
    zlabel(ax, 'Elevation [km]');
    box(ax, 'on');
    grid(ax, 'off');
    set(ax, 'Layer', 'top');
    view(ax, 3);
end

function format_2d_axes(ax, proj, label)
% FORMAT_2D_AXES
% -------------------------------------------------------------------------
    set(ax, 'FontName', 'Verdana', 'FontSize', 12);
    axis(ax, 'equal', 'tight');
    switch lower(proj)
        case 'x'
            xlabel(ax, 'Northing [km]');
            ylabel(ax, 'Elevation [km]');
        case 'y'
            xlabel(ax, 'Easting [km]');
            ylabel(ax, 'Elevation [km]');
        case 'z'
            xlabel(ax, 'Easting [km]');
            ylabel(ax, 'Northing [km]');
        case 'point'
            xlabel(ax, 'Distance along profile [km]');
            ylabel(ax, 'Elevation [km]');
    end

    title(ax, label, 'Interpreter', 'none');
    box(ax, 'on');
    grid(ax, 'off');
    set(ax, 'Layer', 'top');

end

function cvals = draw_plane_3d(ax, P, C, values, P0, n, hexEdges, proj, pval, globalScale)
% DRAW_PLANE_3D
% -------------------------------------------------------------------------
% Draws the polygonal intersection of a plane with all hex cells in 3D.
% Uses hard snapping for axis-aligned planes to reduce cracks.
% -------------------------------------------------------------------------

    nCells = size(C, 1);
    cvals  = [];

    for k = 1:nCells

        vertsIdx = C(k, :);
        verts    = P(vertsIdx, :);

        d = (verts - P0) * n.';

        if all(d > 0) || all(d < 0)
            continue;
        end

        [ptsSlice, valSlice] = intersect_hex_edges( ...
            verts, vertsIdx, values, d, hexEdges, globalScale);

        if size(ptsSlice,1) < 3
            continue;
        end

        % Hard snap for axis-aligned planes
        switch lower(proj)
            case 'x', ptsSlice(:,1) = pval;
            case 'y', ptsSlice(:,2) = pval;
            case 'z', ptsSlice(:,3) = pval;
        end

        if size(ptsSlice,1) < 3
            continue;
        end

        ptsSlice = order_polygon_in_plane(ptsSlice, P0, n);

        cVal = mean(valSlice, 'omitnan');
        cvals(end+1) = cVal; %#ok<AGROW>

        patch(ax, ...
              'XData', ptsSlice(:,1), ...
              'YData', ptsSlice(:,2), ...
              'ZData', ptsSlice(:,3), ...
              'FaceVertexCData', cVal * ones(size(ptsSlice,1),1), ...
              'FaceColor', 'flat', ...
              'EdgeColor', 'none');
    end

end

function cvals = draw_plane_2d(ax, P, C, values, P0, n, hexEdges, proj, pval, plane, globalScale)
% DRAW_PLANE_2D
% -------------------------------------------------------------------------
% Draws the polygonal intersection of a plane with all hex cells projected
% into 2D coordinates. Axis-aligned planes are projected directly to their
% natural coordinate pairs; oblique planes use a local plane basis.
% -------------------------------------------------------------------------

    nCells = size(C, 1);
    cvals  = [];

    for k = 1:nCells

        vertsIdx = C(k, :);
        verts    = P(vertsIdx, :);

        d = (verts - P0) * n.';

        if all(d > 0) || all(d < 0)
            continue;
        end

        [ptsSlice, valSlice] = intersect_hex_edges( ...
            verts, vertsIdx, values, d, hexEdges, globalScale);

        if size(ptsSlice,1) < 3
            continue;
        end

        % Hard snap for axis-aligned planes
        switch lower(proj)
            case 'x', ptsSlice(:,1) = pval;
            case 'y', ptsSlice(:,2) = pval;
            case 'z', ptsSlice(:,3) = pval;
        end

        if size(ptsSlice,1) < 3
            continue;
        end

        % Project into plane coordinates
        switch proj
            case 'x'
                % x = const -> plot (y,z)
                u = ptsSlice(:,2);
                v = ptsSlice(:,3);
        
            case 'y'
                % y = const -> plot (x,z)
                u = ptsSlice(:,1);
                v = ptsSlice(:,3);
        
            case 'z'
                % z = const -> plot (x,y)
                u = ptsSlice(:,1);
                v = ptsSlice(:,2);
        
            case 'point'
                % ---------------------------------------------
                % True profile coordinates:
                %   u = distance along A->B (origin at A)
                %   v = elevation (z)
                % ---------------------------------------------
                A0 = plane.A;
                e  = plane.eProfile;    % unit horizontal direction A->B
                L  = plane.L;
        
                % Vector from A to each intersection point (horizontal component)
                qxy = [ptsSlice(:,1) - A0(1), ...
                       ptsSlice(:,2) - A0(2), ...
                       zeros(size(ptsSlice,1),1)];
        
                % Distance along profile (u = 0 at A)
                u = qxy * e(:);
        
                % Vertical axis
                v = ptsSlice(:,3);
        
                % Optional: discard polygons fully outside [0, L]
                if max(u) < 0 || min(u) > L
                    continue;
                end
        
            otherwise
                % Oblique planes (fallback)
                [e1, e2] = plane_basis(n);
                q = ptsSlice - P0;
                u = q * e1;
                v = q * e2;
        end

        % Order polygon in 2D
        cu = mean(u);
        cv = mean(v);
        ang = atan2(v - cv, u - cu);
        [~, ord] = sort(ang);

        u = u(ord);
        v = v(ord);
        valSlice = valSlice(ord);

        cVal = mean(valSlice, 'omitnan');
        cvals(end+1) = cVal; %#ok<AGROW>

        patch(ax, ...
              'XData', u, ...
              'YData', v, ...
              'FaceVertexCData', cVal * ones(numel(u),1), ...
              'FaceColor', 'flat', ...
              'EdgeColor', 'none');
    end

end

function [ptsSlice, valSlice] = intersect_hex_edges(verts, vertsIdx, values, d, hexEdges, globalScale)
% INTERSECT_HEX_EDGES
% -------------------------------------------------------------------------
% Computes intersection points of the plane with all edges of a hexahedron
% and linearly interpolates scalar values.
%
% A global tolerance is used to reduce numerical cracks across cells.
% -------------------------------------------------------------------------

    ptsSlice = [];
    valSlice = [];

    tol = max(1e-9 * globalScale, 1e-12);

    for e = 1:size(hexEdges,1)

        i = hexEdges(e,1);
        j = hexEdges(e,2);

        di = d(i);
        dj = d(j);

        pi = verts(i,:);
        pj = verts(j,:);

        vi = values(vertsIdx(i));
        vj = values(vertsIdx(j));

        % Both endpoints on plane
        if abs(di) < tol && abs(dj) < tol
            ptsSlice = [ptsSlice; pi]; %#ok<AGROW>
            valSlice = [valSlice; vi]; %#ok<AGROW>

        % First endpoint on plane
        elseif abs(di) < tol
            ptsSlice = [ptsSlice; pi]; %#ok<AGROW>
            valSlice = [valSlice; vi]; %#ok<AGROW>

        % Second endpoint on plane
        elseif abs(dj) < tol
            ptsSlice = [ptsSlice; pj]; %#ok<AGROW>
            valSlice = [valSlice; vj]; %#ok<AGROW>

        % Edge crosses plane
        elseif di * dj < 0
            t = di / (di - dj);
            pInt = pi + t*(pj - pi);
            vInt = vi + t*(vj - vi);

            ptsSlice = [ptsSlice; pInt]; %#ok<AGROW>
            valSlice = [valSlice; vInt]; %#ok<AGROW>
        end
    end

    % Tolerance-aware unique
    if ~isempty(ptsSlice)
        tolP = max(1e-8 * globalScale, 1e-12);
        ptsKey = round(ptsSlice / tolP) * tolP;
        [~, ia] = unique(ptsKey, 'rows', 'stable');
        ptsSlice = ptsSlice(ia, :);
        valSlice = valSlice(ia);
    end

end

function [e1, e2] = plane_basis(n)
% PLANE_BASIS
% -------------------------------------------------------------------------
% Constructs two orthonormal vectors spanning the plane defined by normal n.
% -------------------------------------------------------------------------

    n_unit = n(:) / norm(n);

    tmp = [1; 0; 0];
    if abs(dot(tmp, n_unit)) > 0.99
        tmp = [0; 1; 0];
    end

    e1 = tmp - dot(tmp, n_unit)*n_unit;
    e1 = e1 / norm(e1);

    e2 = cross(n_unit, e1);
    e2 = e2 / norm(e2);

end

function ptsOrdered = order_polygon_in_plane(pts, P0, n)
% ORDER_POLYGON_IN_PLANE
% -------------------------------------------------------------------------
% Orders polygon vertices by angle in local plane coordinates.
% -------------------------------------------------------------------------

    [e1, e2] = plane_basis(n);

    q = pts - P0;
    u = q * e1;
    v = q * e2;

    cu = mean(u);
    cv = mean(v);

    ang = atan2(v - cv, u - cu);
    [~, ord] = sort(ang);

    ptsOrdered = pts(ord, :);

end