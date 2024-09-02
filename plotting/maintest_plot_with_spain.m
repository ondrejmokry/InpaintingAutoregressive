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
methods = ["extrapolation", "janssen", "janssen_hann", "janssen_rect", "aspain", "aspainmod"];

% candidates to plot
orders = [2048, 2048, 1024, 512];
algos = ["arburg", "arburg", "arburg", "arburg"];

% names for the legend
methodnames = [...
    sprintf("extrapolation-based, p = %d, %s", orders(1), algos(1)), ...
    sprintf("Janssen, gap-wise, p = %d, %s", orders(2), algos(2)), ...
    sprintf("Janssen, Hann window, p = %d, %s", orders(3), algos(3)), ...
    sprintf("Janssen, rect. window, p = %d, %s", orders(4), algos(4)), ...
    "A-SPAIN", ...
    "A-SPAIN-MOD" ...
    ];

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
    for m = 1:length(methods)-2 % only AR methods
        tables.(methods(m)) = [tables.(methods(m)); S.tables.(methods(m))];
    end

end
clear a maxit method p S w

%% add SPAIN
S = load(fold + "/results_spain.mat");
tables.aspain = S.tables.aspain;
tables.aspainmod = S.tables.aspainmod;
clear S

for i = 1:length(metrics)

    %% prepare figure
    figure
    colors = colororder;

    %% process
    % dimensions of data
    signals = unique(tables.(methods(1)).signal);
    gaps = unique(tables.(methods(1)).gap);
    data = NaN(length(signals), length(gaps), length(methods));

    % fill the array
    for m = 1:length(methods)
        for s = 1:length(signals)
            for g = 1:length(gaps)
    
                % find the row
                rows = strcmp(tables.(methods(m)).signal, signals(s));
                rows = rows .* (tables.(methods(m)).gap == gaps(g));
                if m <= 4 % only AR methods
                    rows = rows .* (tables.(methods(m)).p == orders(m));
                    rows = rows .* (tables.(methods(m)).method == algos(m));
                end
                row = find(rows);

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
    hold on
    h = gobjects(length(methods), 1);
    for m = 1:length(methods)
        if intervals
            h(m) = fillinterval(gaps, means(:, m), lowers(:, m), uppers(:, m), colors(m, :));
        else
            h(m) = plot(gaps, means(:, m), "Color", colors(m, :));
        end
        % h(m).LineStyle = linestyles(m);
        % h(m).LineWidth = linewidths(m);
    end
    xlim([gaps(1), gaps(end)])
    grid on
    box on
    if i == 1
        legend(h, methodnames, "Location", "northeast")
    else
        legend(h, methodnames, "Location", "southwest")
    end
    title(ftitles(i))
    xlabel("gap length (ms)")
    ylabel(ylabels(i))
    if i > 1
        ylim([-4, 0])
    end
end