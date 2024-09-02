clear
clc
close all

addpath("../utils")

%% settings
bootstrap = true;
intervals = false;
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
ylabels = ["SDR (dB)", "ODG", "ODG"];
ftitles = ["peak SDR", "peak ODG by PEMO-Q", "peak ODG by PEAQ"];

for i = 1:length(metrics)

    fprintf("Metric: %s\n", metrics(i))

    %% prepare figure
    figure
    colors = colororder;
    tls = tiledlayout(2, 5, "TileIndexing", "columnmajor");
    title(tls, ftitles(i))
    
    %% process
    for f = 1:length(filestoload)
    
        %% load data
        fprintf("Loading %s...\n", filestoload(f))
        S = load(fold + "/" + filestoload(f));
        for m = 1:length(methods)
            tables.(methods(m)) = S.tables.(methods(m));
        end
        
        groupname = sprintf("%s, p = %d, w = %d, a = %d", ...
            S.method, S.p, S.w, S.a);
        
        % dimensions of data:
        % number of signals × number of gap lengths × number of methods
        signals = unique(S.tables.(methods(1)).signal);
        gaps = unique(S.tables.(methods(1)).gap);
        data = NaN(length(signals), length(gaps), length(methods));
        for s = 1:length(signals)
            for g = 1:length(gaps)
                for m = 1:length(methods)
    
                    % find the row
                    rows = strcmp(tables.(methods(m)).signal, signals(s));
                    rows = rows .* (tables.(methods(m)).gap == gaps(g));
                    row  = find(rows);
                    if isempty(row)
                        continue
                    end
                    
                    % find maximum
                    data(s, g, m) = max(tables.(methods(m)).(metrics(i)){row});
    
                end
            end
        end
    
        %% interval estimate
        if intervals
            if bootstrap
                % interval estimate using bootstrapping
                fprintf("Computing the bootstraps...\n") %#ok<*UNRCH>
                [means, lowers, uppers] = bootstrap_est(data);
            else
                % interval estimate assuming normality of the data
                fprintf("Computing the stds...\n")
                mult   = tinv(0.975, length(signals)-1);
                means  = squeeze(mean(data, 1, "omitnan"));
                stds   = squeeze(std(data, 0, 1, "omitnan"));
                lowers = means - mult*stds/sqrt(length(signals));
                uppers = means + mult*stds/sqrt(length(signals));
            end
        else
            % no interval estimate
            means = squeeze(mean(data, 1, "omitnan"));
        end
           
        %% plot
        nexttile(tls)
        h = gobjects(length(methods), 1);
        for m = 1:length(methods)
            if intervals
                h(m) = fillinterval(gaps, means(:, m), lowers(:, m), uppers(:, m), colors(m, :));
            else
                h(m) = plot(gaps, means(:, m));
                hold on
            end
            % h(m).LineStyle = linestyles(m);
            % h(m).LineWidth = linewidths(m);
        end
        xlim([gaps(1), gaps(end)])
        grid on
        if i == 1
            legend(h, strrep(methods, "_", " "), "Location", "northeast")
        else
            legend(h, strrep(methods, "_", " "), "Location", "southwest")
        end
        title(groupname)
        xlabel("gap length (ms)")
        ylabel(ylabels(i))
        if i > 1
            ylim([-4, 0])
        end
    end
    
    linkaxes(tls.Children(2:2:end), "xy")

end