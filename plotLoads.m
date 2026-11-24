function plotLoads(model, loads, ax)
%PLOTLOADS  Draw the applied loads on top of the frame drawn by plotModel.
%   plotLoads(model, loads)       draws in the current axes
%   plotLoads(model, loads, ax)   draws in the given axes (e.g. app.FrameAxes)
%
%   red   = tank weight q on the beam
%   blue  = CFD wave drag F_D at the top-left node
%   purple= Morison wave load on the columns (arrow length ~ load size)
%   Call plotModel first, then plotLoads on the same axes.

if nargin < 3 || isempty(ax), ax = gca; end
hold(ax, 'on');
X = model.nodes(:,1); Y = model.nodes(:,2);
W = model.W; H = model.H;
red = [0.85 0.2 0.2]; blue = [0.1 0.35 0.75]; purple = [0.45 0.3 0.75];

% ---- tank weight q: row of down arrows above the beam ----
if loads.q > 0
    hq = 0.35 * H;                                     % arrow length
    xs = linspace(0, W, 11);
    plot(ax, [0 W], (H + hq) * [1 1], '-', 'Color', red, 'LineWidth', 1.2);
    for x = xs
        drawArrow(ax, x, H + hq, x, H, red, 1.2, 0.012*W);
    end
    text(ax, W/2, H + hq, sprintf('q = %.0f kN/m', loads.q/1e3), 'Color', red, ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'bottom', 'FontSize', 10);
end

% ---- CFD vertical force Fz (upward, buoyancy + wave) under the beam ----
if isfield(loads, 'qUp') && loads.qUp > 0
    green = [0.1 0.55 0.35];
    hu = 0.25 * H;
    xc = X(model.topNodes);                           % column positions: keep arrows off them
    xs = linspace(0, W, 23); xs = xs(2:end-1);
    xs = xs(min(abs(xs(:) - xc(:).'), [], 2) > 0.06*W);
    for x = xs
        drawArrow(ax, x, H - hu, x, H - 0.02*H, green, 1, 0.01*W);
    end
    text(ax, mean(xc(1:2)), H - hu, sprintf('F_z = %.0f kN/m up', loads.qUp/1e3), ...
        'Color', green, 'HorizontalAlignment', 'center', 'VerticalAlignment', 'top', 'FontSize', 9);
end

% ---- CFD drag F_D at the top-left node ----
n = model.topNodes(1);
if loads.FD ~= 0
    len = 0.12 * W;
    drawArrow(ax, X(n) - len, Y(n), X(n), Y(n), blue, 2, 0.015*W);
    text(ax, X(n) - len, Y(n), sprintf('F_D = %.0f kN  ', loads.FD/1e3), 'Color', blue, ...
        'HorizontalAlignment', 'right', 'FontSize', 10);
end

% ---- Morison load on the columns ----
wx = loads.wGlobal(:,1);
cols = model.colEls(wx(model.colEls) ~= 0);
if ~isempty(cols)
    wmax = max(abs(wx(cols)));
    lenMax = 0.05 * W;
    for e = cols
        i = model.conn(e,1); j = model.conn(e,2);
        ym = mean(Y([i j])); xm = X(i);
        len = lenMax * wx(e) / wmax;
        drawArrow(ax, xm - abs(len), ym, xm, ym, purple, 1, 0.008*W);
    end
    text(ax, -0.02*W, -0.12*H, sprintf('Morison: %.1f kN total  ', loads.morisonTotal/1e3), ...
        'Color', purple, 'HorizontalAlignment', 'right', 'FontSize', 10);
end

if isfield(loads, 'placeholder') && loads.placeholder
    text(ax, W/2, -0.25*H, 'PLACEHOLDER storm data', 'Color', [0.6 0 0], ...
        'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
end
axis(ax, 'equal'); axis(ax, 'off');
hold(ax, 'off');
end

% ======================================================================
function drawArrow(ax, x0, y0, x1, y1, col, lw, head)
% Straight arrow from (x0,y0) to (x1,y1) with a fixed-size head.
plot(ax, [x0 x1], [y0 y1], '-', 'Color', col, 'LineWidth', lw);
ang = atan2(y1 - y0, x1 - x0);
a1 = ang + pi - pi/7;  a2 = ang + pi + pi/7;
plot(ax, [x1, x1 + head*cos(a1)], [y1, y1 + head*sin(a1)], '-', 'Color', col, 'LineWidth', lw);
plot(ax, [x1, x1 + head*cos(a2)], [y1, y1 + head*sin(a2)], '-', 'Color', col, 'LineWidth', lw);
end