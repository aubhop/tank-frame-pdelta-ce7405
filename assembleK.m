%assembleK
function [K, KE, KG] = assembleK(model, N)
nEl = size(model.conn, 1);
if nargin < 2 || isempty(N)
    N = zeros(nEl, 1);
end
KE = zeros(model.nDOF);
KG = zeros(model.nDOF);

for e = 1:nEl
    sec = model.sections(model.sec(e));
    L = model.L(e);
    T = transformation(model.c(e), model.s(e));
    d = model.dofMap(e,:);
    kLoc = elementStiffness(sec.E, sec.A, sec.I, L, model.type(e));
    KE(d,d) = KE(d,d) + T.' * kLoc * T;
    if N(e) ~= 0
        kgLoc = geometricStiffness(N(e), L, model.type(e));
        KG(d, d) = KG(d, d) + T.' * kgLoc * T;
    end
end
K = KE + KG;
end

function T = transformation(c, s)
lam = [ c s 0
    -s c 0
    0 0 1];
T = blkdiag(lam, lam);
end

function k = elementStiffness(E, A, I, L, type)
a = E*A/L;
k = zeros(6);
k([1 4], [1 4]) = a * [1 -1; -1 1]; %axial terms
if type == 1
    b1 = 12*E*I/L^3;
    b2 = 6*E*I/L^2;
    b3 = 4*E*I/L;
    b4 = 2*E*I/L;
    k([2 3 5 6], [2 3 5 6]) = [b1 b2 -b1 b2; b2 b3 -b2 b4; ...
        -b1 -b2 b1 -b2; b2 b4 -b2 b3];
end
end

function kg = geometricStiffness(N, L, type)
kg = zeros(6);
if type == 1
     kg([2 3 5 6], [2 3 5 6]) = (N/L) * [6/5 L/10 -6/5 L/10
         L/10 2*L^2/15 -L/10 -L^2/30
         -6/5 -L/10 6/5 -L/10
         L/10 -L^2/30 -L/10 2*L^2/15];
else
    kg([2 5], [2 5]) = (N/L) * [1 -1; -1 1];
end
end
