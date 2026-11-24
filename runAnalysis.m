function runAnalysis(app)
%RUNANALYSIS  Everything the app's RUN button does.
%   Called from the app's RUN button callback with one line:
%       runAnalysis(app);
%
%   Components it uses (Component Browser names):
%     Slider, StormDropDown, ContentsDropDown, FillLevelSlider,
%     stOrderCheckBox, ndOrderCheckBox, UIAxes,
%     FrameAxes, DeflectedAxes, ResultsTable, SummaryLabel, StatusLabel
app.StatusLabel.Text = 'Running analysis...';
drawnow;

try
    % ---------------- 1. read the inputs ----------------
    nCols    = round(app.Slider.Value);
    app.Slider.Value = nCols;                    % snap the slider to a whole number
    storm    = app.StormDropDown.Value;
    contents = app.ContentsDropDown.Value;
    fillFrac = app.FillLevelSlider.Value / 100;
    do1 = app.stOrderCheckBox.Value;
    do2 = app.ndOrderCheckBox.Value;
    if ~do1 && ~do2, do1 = true; end             % nothing ticked: show first order

    % ---------------- 2. build the model and loads ----------------
    model = buildFrame(nCols);
    loads = buildLoads(model, storm, contents, fillFrac);

    % ---------------- sections and supports panel ----------------
    sec = model.sections;
    app.SectionsLabel.Text = { ...
        sprintf('Columns: E = %.0f GPa, A = %.3f m², I = %.1e m⁴', sec(1).E/1e9, sec(1).A, sec(1).I), ...
        sprintf('Beam:    E = %.0f GPa, A = %.3f m², I = %.1e m⁴', sec(2).E/1e9, sec(2).A, sec(2).I), ...
        sprintf('Braces:  E = %.0f GPa, A = %.3f m² (axial only)', sec(3).E/1e9, sec(3).A), ...
        'Supports: fixed bases (u = v = θ = 0)'};
    % ---------------- 3. analyses ----------------
    res2 = [];
    lamMax = 1;
    if do2
        res2 = solvePDelta(model, loads);                    % results at lambda = 1 (summary, table)
        lamMax = max(1, 0.8 * res2.lambdaCr);                % plot curves up to 80% of buckling
        curve2 = solvePDelta(model, loads, lamMax, 40);      % extra run just for the curve
        res2.lambda = curve2.lambda; res2.sway = curve2.sway;
    end
    res1 = solveLinear(model, loads, lamMax);                % straight line to the same lambda
    % ---------------- 4. frame view with loads ----------------
    plotModel(model, app.FrameAxes);
    plotLoads(model, loads, app.FrameAxes);
    set(app.FrameAxes, 'DataAspectRatio', [1 0.4 1], ...      % stretch vertical 2.5x
        'XLim', [-0.3 1.08]*model.W, 'YLim', [-0.35 1.9]*model.H);
    title(app.FrameAxes, sprintf('Frame model and loads (%d columns, vertical scale x2.5)', nCols));

    % ---------------- 5. lambda vs. top sway ----------------
    ax = app.UIAxes;
    cla(ax); hold(ax, 'on');
    leg = {};
    if do1
        plot(ax, res1.sway*1e3, res1.lambda, '--', 'Color', [0.1 0.35 0.75], 'LineWidth', 1.5);
        leg{end+1} = '1st order';
    end
    if do2
        plot(ax, res2.sway*1e3, res2.lambda, '-', 'Color', [0.8 0.15 0.15], 'LineWidth', 2);
        leg{end+1} = '2nd order';
    end
    hold(ax, 'off');
    xlabel(ax, 'Top sway (mm)');
    ylabel(ax, 'Load factor \lambda');
    title(ax, '');
    grid(ax, 'on');
    if ~isempty(leg), legend(ax, leg, 'Location', 'southeast'); end

    % ---------------- 6. deflected shapes ----------------
    if do1 && do2
        plotDeflected(model, res1, res2, app.DeflectedAxes);
    elseif do2
        plotDeflected(model, res2, [], app.DeflectedAxes);
    else
        plotDeflected(model, res1, [], app.DeflectedAxes);
    end

    % ---------------- 7. results table: member end forces ----------------
    if do2, r = res2; which = '2nd order'; else, r = res1; which = '1st order'; end
    nEl  = size(model.conn, 1);
    data = [(1:nEl).', model.conn, r.f(:,1)/1e3, r.f(:,2)/1e3, r.f(:,3)/1e3, ...
            r.f(:,4)/1e3, r.f(:,5)/1e3, r.f(:,6)/1e3];
    app.ResultsTable.Data = round(data, 1);
    app.ResultsTable.ColumnName = {'Element', 'Node i', 'Node j', ...
    'Axial force i (kN)', 'Shear i (kN)', 'Bending moment i (kN·m)', ...
    'Axial force j (kN)', 'Shear j (kN)', 'Bending moment j (kN·m)'};
    app.ResultsTable.ColumnFormat = repmat({'shortG'}, 1, 9);

    % ---------------- 8. summary ----------------
    lines = {sprintf('Storm: %s (CFD run #%d)', storm, loads.stormData.simID), ...
             sprintf('Tank q = %.0f kN/m down, F_z = %.0f kN/m up, F_D = %.0f kN', ...
                     loads.q/1e3, loads.qUp/1e3, loads.FD/1e3), ...
             sprintf('Sway, 1st order: %.3f mm', res1.swayTop*1e3)};
    if do2
        AF = res2.swayTop / res1.swayTop;
        lines{end+1} = sprintf('Sway, 2nd order: %.3f mm', res2.swayTop*1e3);
        lines{end+1} = sprintf('Amplification factor AF = %.4f', AF);
        lines{end+1} = sprintf('Critical load factor = %.1f', res2.lambdaCr);
        lines{end+1} = sprintf('NR iterations: %d total', sum(res2.iters));
    end
    lines{end+1} = sprintf('Table shows %s member end forces', which);
    app.SummaryLabel.Text = lines;

    % ---------------- 9. status bar ----------------
    msg = sprintf('%d nodes, %d elements, %d free DOFs. ', ...
        size(model.nodes, 1), nEl, numel(model.freeDOF));
    if do2, msg = [msg res2.message]; else, msg = [msg 'First-order analysis complete.']; end
    if isfield(loads, 'placeholder') && loads.placeholder
        msg = [msg '  WARNING: placeholder storm data.'];
    end
    app.StatusLabel.Text = msg;

catch err
    app.StatusLabel.Text = ['Error: ' err.message];
    uialert(app.UIFigure, err.message, 'Analysis error');
end
end