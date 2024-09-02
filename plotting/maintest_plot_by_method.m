clear
clc
close all

addpath("../utils")

%% settings
bootstrap = false;
intervals = false;
fold = "../results";

filestoload = [ ...
    "results_01", "results_02", ...
    "results_03", "results_04", ...
    "results_05", "results_06", ...
    "results_07", "results_08", ...
    "results_09", "results_10" ...
    ];

if contains(fold, "irmas")
    filestoload = filestoload(1:5);
end

% methods corresponding to fieldnames of the tables variable
methods = ["extrapolation", "janssen", "janssen_hann", ... "janssen_tukey", 
    "janssen_rect"];

% metrics corresponding to variable names of the tables
metrics = ["SDR", "PEMOQ", "PEAQ"];
ylabels = ["SDR (dB)", "ODG", "ODG"];
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

for i = 1:length(metrics)

    fprintf("Metric: %s\n", metrics(i))

    %% prepare figure
    figure
    colors = colororder;
    tls = tiledlayout(1, length(methods));
    title(tls, {ftitles(i), "darker shade: arburg, lighter shade: lpc"})
    
    %% process
    % each method has its own subplot
    for m = 1:length(methods)

        % dimensions of data
        signals = unique(tables.(methods(m)).signal);
        gaps = unique(tables.(methods(m)).gap);
        ps = unique(tables.(methods(m)).p);
        data = NaN(length(signals), length(gaps), length(ps), 2);
        opti = NaN(length(signals), length(gaps), length(ps), 2);
        for s = 1:length(signals)
            for g = 1:length(gaps)
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
                    [data(s, g, p, 1), opti(s, g, p, 1)] = max(tables.(methods(m)).(metrics(i)){burgrow});
                    if lpcrow
                        [data(s, g, p, 2), opti(s, g, p, 2)] = max(tables.(methods(m)).(metrics(i)){lpcrow});
                    end
                    
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
        hold on
        h = gobjects(length(ps), 1);
        for p = 1:length(ps)
            if intervals
                h(p) = fillinterval(gaps, means(:, p, 1), lowers(:, p, 1), uppers(:, p, 1), colors(p, :));
                fillinterval(gaps, means(:, p, 2), lowers(:, p, 2), uppers(:, p, 2), 0.75 + 0.25*colors(p, :))
            else
                h(p) = plot(gaps, means(:, p, 1), "Color", colors(p, :));
                plot(gaps, means(:, p, 2), "Color", 0.6 + 0.4*colors(p, :), "LineStyle", "--")
            end
            % h(p).LineStyle = linestyles(p);
            % h(p).LineWidth = linewidths(p);
        end
        xlim([gaps(1), gaps(end)])
        grid on
        box on
        if i == 1
            legend(h, string(ps), "Location", "northeast")
        else
            legend(h, string(ps), "Location", "southwest")
        end
        title(strrep(methods(m), "_", " "))
        xlabel("gap length (ms)")
        ylabel(ylabels(i))
        if i > 1
            ylim([-4, 0])
        end
    end
    
    linkaxes(tls.Children(2:2:end), "xy")

end