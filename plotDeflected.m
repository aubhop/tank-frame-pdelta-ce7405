function scale = plotDeflected(model, res1, res2, ax, scale)
%PLOTDEFLECTED  Draw the undeformed frame and the first- and second-order deflected shapes.
%
%   plotDeflected(model, res1)              first-order only, new figure
%   plotDeflected(model, res1, res2)        first and second order
%   plotDeflected(model, res1, res2, ax)    draw in given axes (e.g. app.DeflectedAxes)
%   plotDeflected(model, res1, res2, ax, scale)   fixed magnification
%   scale = plotDeflected(...)              returns the magnification used
%
%   Both shapes use the SAME magnification, so they can be compared directly.
%   Frame elements are drawn as smooth curves (cubic beam shape functions).
%
%   gray = undeformed,  blue dashed = first order,  red solid = second order (P-Delta)

if nargin < 3, res2 = []; end
if nargin < 4 || isempty(ax), figure('Color', 'w'); ax = gca; end

% one magnification for both shapes: biggest translation drawn as 8% of the width
uAll = res1.u; if ~isempty(res2), uAll = [uAll, res2.u]; end
trans = uAll(sort([1:3:model.nDOF, 2:3:model.nDOF]), :);
if nargin < 5 || isempty(scale)
    dmax = max(abs(trans(:)));
    if dmax == 0, scale = 1; else, scale = 0.08 * model.W / dmax; end
end

cla(ax); hold(ax, 'on');
gray = [0.7 0.7 0.7]; blue = [0.1 0.35 0.75]; red = [0.8 0.15 0.15];

h0 = drawShape(ax, model, zeros(model.nDOF, 1), 0, gray, '-', 2);
h1 = drawShape(ax, model, res1.u, scale, blue, '--', 1.5);
hs = [h0 h1]; names = {'Undeformed', 'First order'};
if ~isempty(res2)
    h2 = drawShape(ax, model, res2.u, scale, red, '-', 2);
    hs(end+1) = h2; names{end+1} = 'Second order (P-\Delta)';
end

% fixed supports
w = 0.04 * model.W;
for n = model.baseNodes
    plot(ax, model.nodes(n,1) + [-w w], model.nodes(n,2) * [1 1], 'k-', 'LineWidth', 2);
end

axis(ax, 'equal'); axis(ax, 'off');
legend(ax, hs, names, 'Location', 'southoutside', 'Orientation', 'horizontal');
title(ax, sprintf('Deflected shape (magnified x%.0f)', scale));
hold(ax, 'off');
end

% ======================================================================
function h = drawShape(ax, model, u, scale, col, ls, lw)
% Draw every element of the frame displaced by scale*u.
nPts = 15;
h = [];
for e = 1:size(model.conn, 1)
    i = model.conn(e,1);
    c = model.c(e); s = model.s(e); L = model.L(e);
    lam = [c s 0; -s c 0; 0 0 1];
    ul = blkdiag(lam, lam) * u(model.dofMap(e, :));        % local [u1 v1 th1 u2 v2 th2]
    if model.type(e) == 1
        x = linspace(0, L, nPts); xi = x / L;
        ua = ul(1)*(1 - xi) + ul(4)*xi;                    % axial: linear
        N1 = 1 - 3*xi.^2 + 2*xi.^3;  N2 = L*(xi - 2*xi.^2 + xi.^3);
        N3 = 3*xi.^2 - 2*xi.^3;      N4 = L*(-xi.^2 + xi.^3);
        v  = N1*ul(2) + N2*ul(3) + N3*ul(5) + N4*ul(6);     % transverse: cubic
    else
        x = [0 L]; ua = [ul(1) ul(4)]; v = [ul(2) ul(5)];    % truss: straight
    end
    X = model.nodes(i,1) + x*c + scale*(ua*c - v*s);
    Y = model.nodes(i,2) + x*s + scale*(ua*s + v*c);
    lwE = lw; if model.type(e) == 2, lwE = max(0.75, lw/2); end
    hh = plot(ax, X, Y, 'LineStyle', ls, 'Color', col, 'LineWidth', lwE);
    if isempty(h), h = hh; end
end
end