clear
clc
close all

addpath("../utils")

%% settings
displabels = false;
fold = "../results";

filestoload = [ ...
    "results_01", "results_02", ...
    "results_03", "results_04", ...
    "results_05", "results_06", ...
    "results_07", "results_08", ...
    "results_09", "results_10" ...
    ];

% methods corresponding to fieldnames of the tables variable
methods = ["extrapolation", "janssen", "janssen_hann", ... "janssen_tukey",
    "janssen_rect"];

% metrics corresponding to variable names of the tables
metrics = ["SDR", "PEMOQ", "PEAQ"];
axislabels = ["SDR (dB)", "ODG", "ODG"];
ftitles = ["peak SDR", "peak ODG by PEMO-Q", "peak ODG by PEAQ"];

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
    "rejection of the null hypothesis that Burg's algorithm and LPC lead\n" + ...
    "to results with the same median, with the alternate hypothesis that\n" + ...
    "Burg's algorithm leads to results with higher median\n")

for i = 1:length(metrics)

    fprintf("Metric: %s\n", metrics(i))

    % table for p-values
    ps = unique(tables.(methods(1)).p);
    pvals = table('Size', [length(methods), length(ps)], ...
        'VariableTypes', repmat("double", [1, length(ps)]), ...
        'VariableNames', string(ps), 'RowNames', methods);

    %% prepare figure
    figure
    colors = colororder;
    tls = tiledlayout(1, length(methods));
    title(tls, ftitles(i))
    
    %% process
    % each method has its own subplot
    for m = 1:length(methods)

        % dimensions of data
        signals = unique(tables.(methods(m)).signal);
        gaps = unique(tables.(methods(m)).gap);
        ps = unique(tables.(methods(m)).p);
        data = NaN(length(signals), length(gaps), length(ps), 2);


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
                    burgrows = rows .* (tables.(methods(m)).method == "arburg");
                    lpcrows = rows .* (tables.(methods(m)).method == "lpc");
                    burgrow = find(burgrows);
                    lpcrow = find(lpcrows);

                    % find maximum
                    data(s, g, p, 1) = max(tables.(methods(m)).(metrics(i)){burgrow});
                    data(s, g, p, 2) = max(tables.(methods(m)).(metrics(i)){lpcrow});
                    
                end
            end
        end
        
        %% plot
        nexttile(tls)
        hold on
        for p = 1:length(ps)
            x = data(:, :, p, 2); % lpc
            y = data(:, :, p, 1); % arburg
            sc = scatter(x(:), y(:), 18, colors(p, :), "DisplayName", sprintf("p = %d", ps(p)));
            sc.MarkerEdgeAlpha = 0.25;

            if displabels
                % add data labels
                text(x(:), y(:), names(:), "FontSize", 10, "Color", colors(p, :)) %#ok<UNRCH>
            end
        end
        grid on
        box on
        title(strrep(methods(m), "_", " "))
        xlabel(axislabels(i) + ", lpc")
        ylabel(axislabels(i) + ", arburg")
        if i == 1
            xlim([0, 35])
            ylim([0, 35])
        else
            xlim([-4, 0])
            ylim([-4, 0])
        end
        axis square

        %% fill table of p-values
        for p = 1:length(ps)
            burgdata = data(:, :, p, 1);
            lpcdata = data(:, :, p, 2);
    
            pvals{m, p} = signrank(burgdata(:), lpcdata(:), "method", "exact", "tail", "right");
        end
    end
    
    lgd = legend("Orientation", "horizontal");
    lgd.Layout.Tile = "south";
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