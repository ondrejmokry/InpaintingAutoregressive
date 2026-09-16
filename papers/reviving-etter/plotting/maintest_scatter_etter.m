clear
clc
close all

plotDir = fileparts(mfilename("fullpath"));
paperDir = fileparts(plotDir);
repoDir = fileparts(fileparts(paperDir));
addpath(fullfile(repoDir, "utils"))

%% settings
displabels = false;
fold = fullfile(paperDir, "results");

filestoload = [ ...
    "results_01", "results_02", ...
    "results_03", "results_04", ...
    "results_05", "results_06", ...
    "results_07", "results_08", ...
    "results_09", "results_10" ...
    ];

if contains(fold, "longer")
    filestoload = [filestoload, "results_11", "results_12"];
end

% methods corresponding to fieldnames of the tables variable
methods = [ ...
    "extrapolation", ...
    "etter", ...
    "etter_plc"
    ];

estim = "arburg";

% metrics corresponding to variable names of the tables
metrics = ["SDR", "PEMOQ", "PEAQ"];
axislabels = ["SDR (dB)", "ODG", "ODG"];
% ftitles = ["peak SDR", "peak ODG by PEMO-Q", "peak ODG by PEAQ"];
ftitles = ["SDR", "ODG by PEMO-Q", "ODG by PEAQ"];

%% load data
fprintf("Loading %s...\n", filestoload(1))
load(fold + "/" + filestoload(1))
for f = 2:length(filestoload)
    
    fprintf("Loading %s...\n", filestoload(f))
    S = load(fold + "/" + filestoload(f));
    for m = 1:length(methods)
        tables.(methods(m)) = [tables.(methods(m)); S.tables.(methods(m))];
    end

end
clear a maxit method p S w

fprintf("\nWilcoxon signed rank test:\n")
fprintf("For each method and model order, the shown p-value indicates the\n" + ...
    "rejection of the null hypothesis that the two methods lead to\n" + ...
    "results with the same median, with the alternate hypothesis that\n" + ...
    "the 1st method leads to results with higher median\n")

for i = 1:length(metrics)

    fprintf("Metric: %s\n", metrics(i))

    % table for p-values
    ps = unique(tables.(methods(1)).p);
    pvals = table('Size', [4, length(ps)], ...
        'VariableTypes', repmat("double", [1, length(ps)]), ...
        'VariableNames', string(ps), ...
        'RowNames', [ ...
            "Etter versus extrapolation", ...
            "extrapolation versus Etter", ...
            "Etter inpainting versus Etter PLC", ...
            "Etter PLC versus Etter inpainting"]);

    %% prepare figure
    figure
    colors = colororder;
    tls = tiledlayout(1, 2);
    title(tls, ftitles(i))
    
    %% process
    % dimensions of data
    signals = unique(tables.(methods(1)).signal);
    gaps = unique(tables.(methods(1)).gap);
    ps = unique(tables.(methods(1)).p);
    data = NaN(length(signals), length(gaps), length(ps), 3);


    siggap = 0;
    names = strings(length(signals), length(gaps));

    for s = 1:length(signals)
        for g = 1:length(gaps)

            % data label
            siggap = siggap + 1;
            names(s, g) = sprintf("%d/%d", s, g);

            for p = 1:length(ps)

                % find the row
                rows = strcmp(tables.(methods(m)).signal, signals(s));
                rows = rows .* (tables.(methods(m)).gap == gaps(g));
                rows = rows .* (tables.(methods(m)).p == ps(p));
                rows = rows .* (tables.(methods(m)).method == estim);
                row = find(rows);

                % find maximum
                data(s, g, p, 1) = max(tables.(methods(1)).(metrics(i)){row});
                data(s, g, p, 2) = max(tables.(methods(2)).(metrics(i)){row});
                data(s, g, p, 3) = max(tables.(methods(3)).(metrics(i)){row});

            end
        end
    end
    
    for test = 1:2

        switch test
            case 1
                ttl = "Etter versus extrapolation";
                xl = "Etter";
                yl = "extrapolation";
                indx = 1;
                indy = 2;
            case 2
                ttl = "Etter inpainting versus Etter PLC";
                xl = "Etter inpainting";
                yl = "Etter PLC";
                indx = 2;
                indy = 3;
        end

        %% plot
        nexttile(tls)
        hold on
        for p = 1:length(ps)
            x = data(:, :, p, indx);
            y = data(:, :, p, indy);
            sc = scatter(x(:), y(:), 18, colors(p, :), "DisplayName", sprintf("p = %d", ps(p)));
            sc.MarkerEdgeAlpha = 0.5;

            if displabels
                % add data labels
                text(x(:), y(:), names(:), "FontSize", 10, "Color", colors(p, :)) %#ok<UNRCH>
            end
        end
        grid on
        box on
        title(ttl)
        xlabel(axislabels(i) + ", " + xl)
        ylabel(axislabels(i) + ", " + yl)
        if contains(ftitles(i), "SDR")
            xlim([0, 35])
            ylim([0, 35])
        end
        if contains(ftitles(i), "ODG")
            ylim([-4, 0])
            ylim([-4, 0])
        end
        axis square

        %% fill table of p-values
        for p = 1:length(ps)
            xdata = data(:, :, p, indx);
            ydata = data(:, :, p, indy);
    
            pvals{2*test-1, p} = signrank(xdata(:), ydata(:), "method", "exact", "tail", "right");
            pvals{2*test, p} = signrank(ydata(:), xdata(:), "method", "exact", "tail", "right");
        end
    end

    lgd = legend;
    lgd.Layout.Tile = "east";
    % linkaxes(tls.Children(2:end), "xy")

    % add diagonals
    for ax = tls.Children(2:end)'
        newmin = min(ax.XLim(1), ax.YLim(1));
        newmax = max(ax.XLim(2), ax.YLim(2));

        line(ax, [newmin, newmax], [newmin, newmax], "Color", 0.85*[1, 1, 1], "DisplayName", "diagonal")
        set(ax, "XLim", [newmin, newmax])
        set(ax, "YLim", [newmin, newmax])
    end

    % display p-values
    disp(pvals)

end