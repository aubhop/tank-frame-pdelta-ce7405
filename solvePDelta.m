function res = solvePDelta(model, loads, lambdaMax, nSteps, tol, maxIter)
if nargin < 3 || isempty(lambdaMax)
    lambdaMax = 1;

end
if nargin < 4 || isempty(nSteps)
    nSteps = 10;

end
if nargin < 5 || isempty(tol)
    tol = 1e-6;

end
if nargin < 6 || isempty(maxIter)
    maxIter = 30;
end

f = model.freeDOF;
s = model.fixedDOF;
nEl = size(model.conn, 1);
F = loads.F;
topDOF = 3*model.topNodes(1) -2;

% load factors: equal steps
lamList = linspace(0, lambdaMax, nSteps + 1);
if lambdaMax >= 1, lamList = unique([lamList, 1]);
end
lamList = lamList(2:end);

[~, KE] = assembleK(model);
u = zeros(model.nDOF, 1);
lamOut = 0;
swaOut = 0;
iters = zeros(1, 0);
converged = true;
message = '';
saved = [];

for lam = lamList
    ok = false;
    for it = 1:maxIter
        N = axialForces(model, loads, u, lam);
        [~, ~, KG] = assembleK(model, N);
        KT = KE + KG;
        R = lam*F(f) - KT(f, f) * u(f);
        if norm(R) <= tol * norm(lam*F(f))
            ok = true; break;
        end
        [Rc, p] = chol(KT(f, f));
        if p > 0
            break;
        end
         u(f) = u(f) + Rc \ (Rc.' \ R);
    end
    if ~ ok
        converged = false;
        message = sprintf('Unstable or not converged at lambda %.3f (last converged lambda = %.3f).', lam, lamOut(end));
        break;
    end
    lamOut(end+1) = lam;
    swaOut(end+1) = u(topDOF);
    iters(end+1) = it;
    if abs(lam - 1) < 1e-12
        saved = struct('u', u, 'KT', KT, 'N', N, 'lam', lam);
    end
end
if isempty(saved)
    %converged state
    N = axialForces(model, loads, u, lamOut(end));
     [~, ~, KG] = assembleK(model, N);
     saved = struct('u', u, 'KT', KE + KG, 'N', N, 'lam', lamOut(end));
end
if converged
    message = sprintf('Converged: %d load steps, %d iterations in total.', numel(iters), sum(iters));
end

%results 

uS = saved.u;
lamS = saved.lam;
Rr = zeros(model.nDOF, 1);
Rr(s) = saved.KT(s, f) * uS(f) - lamS * F(s);
fEnd = zeros(nEl, 6);
for e = 1:nEl
    fEnd(e,:) = elementEndForces(model, loads, uS, lamS, e, saved.N(e));
end

res.u = uS;
res.R = Rr;
res.f = fEnd;
res.N = -fEnd(:,1);
res.swayTop = uS(topDOF);
res.lambda = lamOut;
res.sway = swaOut;
res.iters = iters;
res.converged = converged;
res.message = message;
res.lambdaCr = criticalLoadFactor(model, loads, KE);
end

function N = axialForces(model, loads, u, lam)
nEl = size(model.conn, 1);
N = zeros(nEl, 1);
for e = 1:nEl
    fe = elementEndForces(model, loads, u, lam, e, 0);
    N(e) = -fe(1);
end
end
function fe = elementEndForces(model, loads, u, lam, e, Ne);
sec = model.sections(model.sec(e));
c = model.c(e);
s = model.s(e);
L = model.L(e);
lamT = [c s 0; -s c 0; 0 0 1];
T = blkdiag(lamT, lamT);
k = zeros(6);
k([1 4], [1 4]) = (sec.E*sec.A/L) * [1 -1; -1 1];
if model.type(e) == 1
    E = sec.E; I = sec.I;
    b1 = 12*E*I/L^3; b2 = 6*E*I/L^2; b3 = 4*E*I/L; b4 = 2*E*I/L;
    k([2 3 5 6], [2 3 5 6]) = [b1 b2 -b1 b2; b2 b3 -b2 b4; -b1 -b2 b1 -b2; b2 b4 -b2 b3];
    if Ne ~= 0
        k([2 3 5 6], [2 3 5 6]) = k([2 3 5 6], [2 3 5 6]) + (Ne/L) * ...
            [6/5 L/10 -6/5 L/10; L/10 2*L^2/15 -L/10 -L^2/30; -6/5 -L/10 6/5 -L/10; L/10 -L^2/30 -L/10 2*L^2/15];
    end
elseif Ne ~= 0
    k([2 5], [2 5]) = (Ne/L) * [1 -1; -1 1];
end
fe = (k * T * u(model.dofMap(e, :))).' - lam * loads.fEqLocal(e, :);
end

%linear buckling estimate
function lamCr = criticalLoadFactor(model, loads, KE)
f = model.freeDOF;
u1 = zeros(model.nDOF, 1);
u1(f) = KE(f, f) \ loads.F(f);
N1 = axialForces(model, loads, u1, 1);
[~, ~, KG1] = assembleK(model, N1);
if ~any(N1 < 0), lamCr = Inf; return; end
ev = eig(KE(f, f), -KG1(f, f));
ev = real(ev(isfinite(ev) & abs(imag(ev)) < 1e-9 & real(ev) > 0));
if isempty(ev), lamCr = Inf; else, lamCr = min(ev); end
end