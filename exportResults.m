function exportResults(app)
fig = ancestor(app.ResultsTable, 'figure');
% EXPORTRESULTS  Save the Results Table and the Summary to an Excel file.
% Called from File > Export Results in CE7405_FinalProj2.

    if isempty(app.ResultsTable.Data)
        uialert(app.TankFramePDeltaAnalyzerUIFigure, 'Press RUN first, then export.', 'No results');
        return
    end

    [file, folder] = uiputfile('*.xlsx', 'Save results as', 'TankFrame_results.xlsx');
    figure(app.TankFramePDeltaAnalyzerUIFigure);                 % bring the app back to the front
    if isequal(file, 0), return; end      % user pressed Cancel
    fname = fullfile(folder, file);

    try
        % Sheet 1: member end forces (what the Results Table shows)
        data = app.ResultsTable.Data;
        if istable(data)
            T = data;
        else
            names = cellstr(app.ResultsTable.ColumnName);
            T = array2table(data, 'VariableNames', names);
        end
        writetable(T, fname, 'Sheet', 'Member forces');

        % Sheet 2: summary text (storm, loads, sway, AF, lambda_cr)
        summ = app.SummaryLabel.Text;
        summ = cellstr(summ);
        writecell(summ(:), fname, 'Sheet', 'Summary');

        app.StatusLabel.Text = ['Results saved to ' file];
    catch err
        uialert(app.TankFramePDeltaAnalyzerUIFigure, err.message, 'Export failed');
    end
end
