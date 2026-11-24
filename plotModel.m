%PLOT MODEL

function plotModel(model, ax)

if nargin < 2 || isempty(ax)
    figure('Color', 'w'); 
    ax = gca;
end
cla(ax);
hold(ax, 'on');
axis(ax, 'equal');
axis(ax, 'off');

X = model.nodes(:,1);
Y = model.nodes(:,2);

for e = 1:size(model.conn, 1)
    i = model.conn(e,1);
    j = model.conn(e, 2);
    if model.type(e) == 2
        plot(ax, X([i j]), Y([i j]), 'k-', 'LineWidth', 1);
    else
        plot(ax, X([i j]), Y([i j]), 'k-', 'LineWidth', 3);
    end 
    t = 0.5;
    if model.type(e) == 2
        t = 1/3;
    end
    xm = X(i) + t*(X(j)-X(i));
    ym = Y(i) + t*(Y(j)-Y(i));
    d = 0.02*model.W;
    %labelling element numbers (blue)
    text(ax, xm - d*model.s(e), ym + d*model.c(e), sprintf('%d', e), ...
        'Color', [0.1 0.3 0.8], 'FontSize', 10, 'HorizontalAlignment', 'center');
end
plot(ax, X, Y, 'k.');

%node loop, going to label node numbers
for n = 1:size(model.nodes, 1)
    text(ax, X(n), Y(n), sprintf(' %d', n), 'FontSize', 10)
end

% FIXED SUPPORTS

w = 0.04 * model.W;
for n = model.baseNodes
    plot(ax, X(n) + [w -w], Y(n) * [1 1], 'k-', 'LineWidth', 2);
end
title(ax, sprintf('%d columns, %d bays, %d nodes, %d elements, %d free DOFs', ...
    model.nCols, model.nCols -1, size(model.nodes, 1), size(model.conn, 1), numel(model.freeDOF)));
hold( ax, 'off');
end