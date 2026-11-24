function loads = buildLoads(model, storm, contents, fillFrac, opts)
%build the load vector for frame

% non changing parameters
% tank radius: 10 m

%legend: H = wave height, d = initial wave depth, T = period, Fx = peak horizontal force on tank



g = 9.81; %grav const

if nargin < 5
    opts = struct();
end

placeholder = false;
stormTable = struct(...
    'name',  {'Moderate', 'Extreme', 'Very Extreme'},...
    'simID', {57, 24, 40},...
    'Fx',    {17.30e3, 98.00e3, 305.19e3},...
    'Fz',    {176.31e3, 565.38e3, 874.39e3},...
    'H',     {1.480, 3.745, 6.808},...
    'T',     {6, 6, 10},...
    'd',     {5, 8, 10});

%------------fixed values---------
r = 10;
Dleg = 1.0;


% ---------------options-----------
s = getOpt(opts, 's', 5); %frame spacing
useMor = getOpt(opts, 'morison', true);
useFz = getOpt(opts, 'buoyancy', true); %CFD upward force on/off
CD = getOpt(opts, 'CD', 1.0); %morison drag
CM  = getOpt(opts, 'CM', 2.0);
rhoW = getOpt(opts, 'rhoW', 1025); %density

if fillFrac < 0 || fillFrac > 1
    error('buildLoads:fill', 'fillFrace must be between 0 and 1.');
end

% ----------storm case-------------
iS = find(strcmpi({stormTable.name}, storm), 1);
if isempty(iS)
    error('buildLoads:storm', 'Unknown storm case "%s".', storm);
end
sd = stormTable(iS);

% --------contents denisty--------
switch lower(contents)
    case 'water'
        rho = 1000;
    case 'liquid hydrogen'
        rho = 70.8;
    otherwise
        error('buildLoads:contents', 'Unknown contents "%s".', contents);
end

nEl = size(model.conn, 1);
wGlobal = zeros(nEl, 2);

%------------tank weight---------
wContents = fillFrac*rho*g*pi*r^2; % N/m of tank length
Wtank = wContents*s; % total on one fram (N)
q = Wtank / model.W; %spread over the beam (N/m)
wGlobal(model.beamEls, 2) = wGlobal(model.beamEls, 2) - q;

%------------CFD vertical force on tank (buoyancy + wave)---------
Fup = 0;
qUp = 0;
if useFz
    Fup = sd.Fz * s; % N on one frame (+ upward)
    qUp = Fup / model.W; % N/m on the beam
    wGlobal(model.beamEls, 2) = wGlobal(model.beamEls, 2) + qUp;
end

%----------------CFD Wave drag at the top left node ------------------
FD = sd.Fx * s;
Fnodal = zeros(model.nDOF, 1);
nTop = model.topNodes(1);
Fnodal(3*nTop - 2) = Fnodal(3*nTop - 2) + FD;

%--------------Morison wave load on the columns--------------------
phase = NaN;
morTotal = 0;
if useMor
    [wCol, phase] = morisonColumns(model, sd, Dleg, CD, CM, rhoW, g);
    wGlobal(:,1)=wGlobal(:,1)+wCol;
    morTotal = sum(wCol .* model.L);
end

%-----------equivalent nodal loads------------
fEqLocal = zeros(nEl, 6);
Feq = zeros(model.nDOF, 1);
for e = 1:nEl
    if all(wGlobal(e, :) == 0) || model.type(e) ==2, continue; 
    end
    c = model.c(e);
    sn = model.s(e);
    L = model.L(e);
    wx = c*wGlobal(e,1) + sn * wGlobal(e,2);
    wy = -sn * wGlobal(e,1) + c *wGlobal(e,2);
    f = [wx*L/2; wy*L/2; wy*L^2/12; wx*L/2; wy*L/2; -wy*L^2/12];
    fEqLocal(e,:) = f.';
    T = transformation(c, sn);
    d = model.dofMap(e, :);
    Feq(d) = Feq(d) + T.' * f;
end

%----------output------
loads.F = Fnodal + Feq;
loads.Fnodal = Fnodal;
loads.Feq = Feq;
loads.wGlobal = wGlobal;
loads.fEqLocal = fEqLocal;
loads.q = q;
loads.FD = FD;
loads.Fup = Fup;
loads.qUp = qUp;
loads.Wtank = Wtank;
loads.rho = rho;
loads.stormData = sd;
loads.morisonPhase = phase;
loads.morisonTotal = morTotal;
loads.placeholder = placeholder;
end

function [wCol, phaseMax] = morisonColumns(model, sd, Dleg, CD, CM, rhoW, g)
nEl = size(model.conn, 1);
wCol = zeros(nEl, 1);
cols = model.colEls;
if isempty(cols) 
    phaseMax = NaN;
    return
end
w = 2*pi /sd.T;
k = fzero(@(k) w^2 - g*k*tanh(k*sd.d), w^2/g); %dispersion relation
yMid = 0.5 * (model.nodes(model.conn(cols,1), 2) + model.nodes(model.conn(cols, 2), 2));
A = cosh( k*yMid) / sinh(k*sd.d); %depth decay
crest = sd.d + sd.H/2;
wet = yMid <= crest;

phases = linspace(0, 2*pi, 361);
best = -inf; 
phaseMax = 0;
fBest = zeros(numel(cols), 1);
for th = phases
    u = (pi*sd.H/sd.T)*A*cos(th);
    dudt = (2*pi^2*sd.H/sd.T^2)*A*sin(th);
    f = 0.5*rhoW*CD*Dleg*u.*abs(u) + rhoW*CM*(pi*Dleg^2/4)*dudt;
    f(~wet) = 0;
    total = sum(f .* model.L(cols));
    if total > best
        best = total;
        phaseMax = th;
        fBest = f;
    end
end
wCol(cols) = fBest;
end

function T = transformation(c, s)
lam = [c s 0; -s c 0; 0 0 1];
T = blkdiag(lam, lam);
end

function v = getOpt(opts, field, default)
if isfield(opts, field) && ~isempty(opts.(field))
    v = opts.(field);
else
    v = default;
end
end
