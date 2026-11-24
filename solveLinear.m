function res = solveLinear(model, loads, lambdaMax, nSteps)
if nargin < 3 || isempty(lambdaMax)
    lambdaMax = 1; 
end
if nargin < 4 || isempty(nSteps)
    nSteps = 10;
end
K = assembleK(model); %elastic stiffness only
f = model.freeDOF;
F = loads.F;
s = model.fixedDOF;


%solve free dofs
u = zeros(model.nDOF, 1);
u(f) = K(f, f)\F(f);

% support reactions
R = zeros(model.nDOF, 1);
R(s) = K(s, f) * u(f) - F(s);

% member end forces
[fEnd, N] = memberForces(model, loads, u, zeros(size(model.conn, 1), 1));

%lambda vs sway
topDOF = 3*model.topNodes(1) - 2;
lambda = linspace(0, lambdaMax, nSteps + 1);

res.u = u;
res.R = R;
res.f = fEnd;
res.N = N;
res.lambda = lambda;
res.sway = lambda * u(topDOF);
res.swayTop = u(topDOF);
res.iters = 1;
end

function [fEnd, N] = memberForces(model, loads, u, Nprev)

nEl = size(model.conn, 1);
fEnd = zeros(nEl, 6);
for e = 1:nEl
    sec = model.sections(model.sec(e));
    c = model.c(e); s = model.s(e); L = model.L(e);
    lam = [ c s 0; -s c 0; 0 0 1];
    T = blkdiag(lam, lam);
    kLoc = localK(sec.E, sec.A, sec.I, L, model.type(e), Nprev(e));
    ue = u(model.dofMap(e,:));
    fEnd(e, :) = (kLoc * T * ue).' - loads.fEqLocal(e, :);
end
N = -fEnd(:,1);
end

function k = localK(E, A, I, L, type, N)
k = zeros(6);
k([1 4], [1 4]) = (E*A/L) * [1 -1; -1 1];
if type == 1
    b1 = 12*E*I/L^3; 
    b2 = 6*E*I/L^2;
    b3 = 4*E*I/L;
    b4 = 2*E*I/L;
    k([2 3 5 6], [2 3 5 6]) = [ b1 b2 -b1 b2;
                                   b2 b3 -b2 b4;
                                  -b1 -b2 b1 -b2;
                                   b2 b4 -b2 b3];
    if N ~= 0
         k([2 3 5 6], [2 3 5 6]) =  k([2 3 5 6], [2 3 5 6]) + (N/L) * ...
             [6/5 L/10 -6/5 L/10; L/10 2*L^2/15 -L/10 -L^2/30; -6/5 -L/10 6/5 -L/10; ...
             L/10 -L^2/30 -L/10 2*L^2/15];
    end
elseif N ~= 0
    k([2 5], [2 5]) = (N/L) * [1 -1; -1 1];
end
end