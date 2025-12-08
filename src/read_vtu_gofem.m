function M = read_vtu_gofem(filename)
% READ_VTU_GOFEM
% -------------------------------------------------------------------------
% Reads a GoFEM VTU (VTK XML Unstructured Grid) file saved in ASCII format 
% and extracts basic mesh information for use in MATLAB.
%
% This simple reader assumes:
%   - Single <UnstructuredGrid><Piece> block
%   - <Points> has a single DataArray with 3 components (x,y,z)
%   - <Cells> has DataArray with Name="connectivity", "offsets", "types"
%   - <PointData> contains a DataArray Name="conductivity" (scalar per node)
%
% OUTPUT
%   M : struct with fields
%       .points      node coordinates (x,y,z) (km)      [Npts x 3]   
%       .cells       connectivity (1-based, NaN padded) [Ncells x maxNodes] 
%       .cellTypes   VTK cell types (integer)           [Ncells x 1] 
%       .pointData   scalar field at nodes
%           .conductivity log10(S/m)                    [Npts x 1] 
%           .resistivity  log10(Ohm·m)                  [Npts x 1]
% -------------------------------------------------------------------------
% UPDATE LOG
%   version         description
%   1.0             First release
% -------------------------------------------------------------------------
%   Author:         César Daniel Castro
%   Affiliation:    Institute of Geophysics, Czech Academy of Science
%   First release:  12-Dec-2025
%   Last update:    12-Dec-2025
% -------------------------------------------------------------------------

    % Read entire file as char
    txt = fileread(filename);

    % POINTS (coordinates)
    % ---------------------------------------------------------------------
    pointsBlock = extractBlock(txt, '<Points>', '</Points>');
    % Inside <Points> we assume a single DataArray
    pointsData  = extractDataArray(pointsBlock);
    coords      = sscanf(pointsData, '%f');      % flat vector
    if mod(numel(coords), 3) ~= 0
        error('Number of coordinate values is not multiple of 3.');
    end
    M.points = reshape(coords, 3, []).';         % [Npts x 3]
    % Coordinate transform
    % The VTU file uses:            We want:
    %   x = North-South                 x = East-West
    %   y = East-West                   y = North-South
    %   z must be inverted (positive <-> negative).
    M.points = [M.points(:,2), M.points(:,1), -M.points(:,3)]./1000;

    % CELLS (connectivity, offsets, types)
    % ---------------------------------------------------------------------
    cellsBlock = extractBlock(txt, '<Cells>', '</Cells>');

    connData = extractDataArrayByName(cellsBlock, 'connectivity');
    if isempty(connData)
        error('Could not find DataArray Name="connectivity" in <Cells>.');
    end
    conn = sscanf(connData, '%d');              % 0-based indices (VTK)

    offsData = extractDataArrayByName(cellsBlock, 'offsets');
    if isempty(offsData)
        error('Could not find DataArray Name="offsets" in <Cells>.');
    end
    offsets = sscanf(offsData, '%d');           % cumulative counts

    typesData = extractDataArrayByName(cellsBlock, 'types');
    if isempty(typesData)
        error('Could not find DataArray Name="types" in <Cells>.');
    end
    M.cellTypes = sscanf(typesData, '%d');

    % Build padded connectivity matrix [Ncells x maxNodesPerCell]
    nCells     = numel(offsets);
    nodeCounts = diff([0; offsets]);            % nodes per cell
    maxNodes   = max(nodeCounts);

    cells = NaN(nCells, maxNodes);
    startIdx = 1;
    for ic = 1:nCells
        endIdx = offsets(ic);
        idx    = conn(startIdx:endIdx);         % 0-based indices
        startIdx = endIdx + 1;

        % Convert to 1-based for MATLAB
        idx = double(idx) + 1;
        cells(ic, 1:numel(idx)) = idx(:).';
    end
    M.cells = cells;

    % POINT DATA (conductivity)
    % ---------------------------------------------------------------------
    pointDataBlock = extractBlock(txt, '<PointData', '</PointData>');
    condData = extractDataArrayByName(pointDataBlock, 'conductivity');
    if isempty(condData)
        warning('No PointData field named "conductivity" found.');
        M.pointData.conductivity = [];
        M.pointData.resistivity  = [];
    else
        % Conductivity (as stored in the VTU) – assumed S/m
        sigma = sscanf(condData, '%f');
        M.pointData.conductivity = log10(sigma);

        % Resistivity = 1 / sigma  (Ohm·m)
        % Basic safeguard to avoid division by zero / negative values
        rho = NaN(size(sigma));
        mask = sigma > 0;                % only positive conductivities
        rho(mask) = 1 ./ sigma(mask);    % ρ = 1/σ
        M.pointData.resistivity = log10(rho);
    end

end

function block = extractBlock(txt, startTag, endTag)
% EXTRACTBLOCK
% -------------------------------------------------------------------------
% Extract first block between <startTag> and </endTag>
% -------------------------------------------------------------------------

    % Use non-greedy match to capture smallest block
    pattern = [regexptranslate('escape', startTag), '.*?', ...
               regexptranslate('escape', endTag)];
    token = regexp(txt, pattern, 'match', 'once');
    if isempty(token)
        error('Could not find block %s ... %s in VTU file.', startTag, endTag);
    end
    block = token;

end

function dataStr = extractDataArray(block)
% EXTRACTDATAARRAY
% -------------------------------------------------------------------------
% Extract content between <DataArray ...> and </DataArray>
% -------------------------------------------------------------------------

    pattern = '<DataArray[^>]*>(.*?)</DataArray>';
    token = regexp(block, pattern, 'tokens', 'once');
    if isempty(token)
        error('No <DataArray> found in the given block.');
    end
    dataStr = token{1};

end

function dataStr = extractDataArrayByName(block, name)
% EXTRACTDATAARRAYBYNAME
% -------------------------------------------------------------------------
% Extract DataArray with a given Name="...", returns empty [] if not found
% -------------------------------------------------------------------------

    pattern = ['<DataArray[^>]*Name="', name, '"[^>]*>(.*?)</DataArray>'];
    token = regexp(block, pattern, 'tokens', 'once');
    if isempty(token)
        dataStr = [];
    else
        dataStr = token{1};
    end

end