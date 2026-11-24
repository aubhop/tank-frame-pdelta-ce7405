%frame
function model = buildFrame(nCols, opts)
%--------- inputs and defualts ---------------------------------
if nargin < 2, opts = struct();
end
if nCols < 2 || nCols ~= round(nCols)
    error('buildFrame:nCols', 'Number of columns must be a whole number >= 2.');
end

W  = getOpt(opts, 'W', 20); % total width
H = getOpt(opts, 'H', 3); % 3 meter column height NEEDS TO BE CHANGED *********
hBrace = getOpt(opts, 'hBrace', H/2);
nElCol = getOpt(opts, 'nElCol', 4);
nElBeam = getOpt(opts, 'nElBeam', 1);
sections = getOpt(opts, 'sections', defaultSections());
if hBrace <= 0 || hBrace > H
    error('buildFrame:hBrace', 'Brace hieght must be between 0 and H.');
end

nBays = nCols - 1;
xCol = linspace(0, W, nCols); %column x position (equal bays)

% -------NODES----------------------
yCol = unique([linspace(0, H, nElCol + 1), hBrace]);
nPerCol = numel(yCol);
nodes = zeros(nCols*nPerCol, 2);
colNode = zeros(nCols, nPerCol);
n = 0;
for c = 1:nCols
    for k = 1:nPerCol
        n = n + 1;
        nodes(n, :) = [xCol(c), yCol(k)];
        colNode(c, k) = n;
    end
end
baseNodes = colNode(:, 1).';
topNodes = colNode(:, end).';
kBrace = find(abs(yCol - hBrace) < 1e-9, 1);
braceNodes = colNode(:, kBrace).';

beamNodeRow = topNodes(1);
for b = 1:nBays
    xs = linspace(xCol(b), xCol(b+1), nElBeam +1);
    for k = 2:nElBeam
        n = n + 1;
        nodes(n, :) = [xs(k), H];
        beamNodeRow(end+1) = n;
    end
    beamNodeRow(end+1) = topNodes(b+1);
end

% -------ELEMENTS-------------------
conn = zeros(0, 2);
type = zeros(0, 1);
sec = zeros(0,1);
group = {};
colOfEl = zeros(0,1);
iCol = findSection(sections, 'COL');
iBeam = findSection(sections, 'BEAM');
iBrace = findSection(sections, 'BRACE');
for c = 1:nCols
    for k = 1:nPerCol-1
        conn(end+1,:) = [colNode(c,k), colNode(c,k+1)];
        type(end+1,1) = 1; sec(end+1,1) = iCol;
        group{end+1,1} = 'COL'; colOfEl(end+1,1) = c;
    end
end

for k = 1:numel(beamNodeRow)-1
    conn(end+1,:) = [beamNodeRow(k), beamNodeRow(k+1)];
    type(end+1,1) = 1;
    sec(end+1, 1) = iBeam;
    group{end+1,1} = 'beam';
    colOfEl(end+1,1) = 0;
end

for b = 1:nBays
    conn(end+1, :) = [baseNodes(b), braceNodes(b+1)];
    conn(end+1, :) = [baseNodes(b+1), braceNodes(b)];
    for r = 1:2
        type(end+1,1) = 2;
        sec(end+1,1) = iBrace;
        group{end+1,1} = 'brace';
        colOfEl(end+1,1) = 0;
    end
end

% --------Geometry of elements ---------
nEl = size(conn, 1);
dx = nodes(conn(:, 2), 1) - nodes(conn(:, 1), 1);
dy = nodes(conn(:, 2), 2) - nodes(conn(:, 1), 2);
L  = hypot(dx, dy);
cs = dx ./ L;
sn = dy ./ L;

% ------dof--------
nNode = size(nodes, 1);
nDOF  = 3 * nNode;
dofMap = zeros(nEl, 6);
for e = 1:nEl
    i = conn(e, 1); j = conn(e, 2);
    dofMap(e, :) = [3*i-2, 3*i-1, 3*i, 3*j-2, 3*j-1, 3*j];
end
fixedDOF = sort(reshape([3*baseNodes-2; 3*baseNodes-1; 3*baseNodes], 1, []));
freeDOF  = setdiff(1:nDOF, fixedDOF);

%-----output-------
model.nodes     = nodes;
model.conn      = conn;
model.type      = type;
model.group     = group;
model.sec       = sec;
model.sections  = sections;
model.L         = L;
model.c         = cs;
model.s         = sn;
model.dofMap    = dofMap;
model.nDOF      = nDOF;
model.fixedDOF  = fixedDOF;
model.freeDOF   = freeDOF;
model.baseNodes = baseNodes;
model.topNodes  = topNodes;
model.beamEls   = find(strcmp(group, 'beam')).';
model.colEls    = find(strcmp(group, 'COL')).';
model.braceEls  = find(strcmp(group, 'brace')).';
model.colOfEl   = colOfEl;
model.W = W; model.H = H; model.hBrace = hBrace; model.nCols = nCols;

end

function s = defaultSections()
% Placeholder steel sections -- replace A and I with your chosen members *** 
E = 200e9;                                       % Pa
s(1) = struct('name', 'COL',   'E', E, 'A', 0.05, 'I', 2.0e-3);
s(2) = struct('name', 'BEAM',  'E', E, 'A', 0.05, 'I', 2.0e-3);
s(3) = struct('name', 'BRACE', 'E', E, 'A', 0.01, 'I', 0);      % truss: I unused
end
 
function idx = findSection(sections, name)
idx = find(strcmpi({sections.name}, name), 1);
if isempty(idx)
    error('buildFrame:section', 'Section "%s" not found in opts.sections.', name);
end
end
 
function v = getOpt(opts, field, default)
if isfield(opts, field) && ~isempty(opts.(field))
    v = opts.(field);
else
    v = default;
end
end