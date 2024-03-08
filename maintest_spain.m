clear
clc
close all

addpath(genpath("utils"))
addpath(genpath("references"))

% future file name
filename = "results_spain";

%% load signals and masks
load("gaps_table.mat")
gaplengths = 10:10:80; % this corresponds to what is saved in gaps_table

%% set params
w = 4096;
a = 1024;
M = 4096;

% SPAIN
SPA.algorithm = 'aspain'; % algorithm
SPA.w = w;                % window length
SPA.a = a;                % window shift
SPA.wtype = 'hann';       % window shape
SPA.F = frame('dft');     % frequency transform
SPA.F.redundancy = 2;     % frequency transform redundancy
SPA.F.frana = @(insig)dft([insig; zeros(length(insig)*(SPA.F.redundancy-1),1)]);
SPA.F.frsyn = @(insig)postpad(idft(insig),length(insig)/SPA.F.redundancy);
SPA.s = 1;                % increment of k
SPA.r = 1;                % every r-th iteration increment k by s   
SPA.epsilon = 0.01;       % stopping criterion of termination function
SPA.maxit = ceil(floor(SPA.w*SPA.F.redundancy/2+1)*SPA.r/SPA.s); % maximum number of iterations
SPA.store_snr = false;    % save SNR during iterations?
SPA.store_obj = false;    % save objective during iterations?

% SPAIN modified
SPM.algorithm = 'aspain_mod'; % formulation
SPM.w = w;                % window length
SPM.a = a;                % window shift
SPM.wtype = 'hann';       % window shape
SPM.M = M;                % number of frequency channels
SPM.gwindow = gabwin(SPM.wtype, SPM.a, SPM.M);
SPM.gwindow = normalize(SPM.gwindow, 'peak'); % peak-normalization of the analysis window
SPM.gdual = gabdual(SPM.gwindow, SPM.a, SPM.M); 
SPM.s = 1;                % increment of k
SPM.r = 1;                % every r-th iteration increment k by s   
SPM.epsilon = 0.01;       % stopping criterion of termination function
SPM.maxit = ceil(floor(SPM.w/2+1)*SPM.r/SPM.s)*2; % maximum number of iterations
SPM.store_snr = false;    % save SNR during iterations?
SPM.store_obj = false;    % save objective during iterations?

%% prepare table for saving everything
methods = ["aspain", "aspainmod"];
types = ["string", "double", "double", "cell", "cell", ...
    "cell", "double", "cell", "cell", "cell"];
names = ["signal", "fs", "gap", "clean", "mask", ...
    "restored", "time", "SDR", "PEMOQ", "PEAQ"];
units = ["", "Hz", "ms", "", "", "", "s", "dB", "", ""];
selection = [1, 2, 3, 7:10]; % all except clean, mask and restored signals
rows = height(gaps_table) * length(gaplengths);
for m = 1:length(methods)
    sigtables.(methods(m)) = table('Size', [length(gaplengths), length(types)], ...
        'VariableTypes', types, 'VariableNames', names);
    sigtables.(methods(m)).Properties.VariableUnits = units;
    tables.(methods(m)) = table('Size', [rows, length(selection)], ...
        'VariableTypes', types(selection), 'VariableNames', names(selection));
    tables.(methods(m)).Properties.VariableUnits = units(selection);
end

% if the experiment has been stopped, load the data
if isfile(filename + ".mat")
    load(filename + ".mat")
end
skipped = 0;

% save time for global timer
cinit = datetime("now");

% start processing
row = 0;
for i = 1:height(gaps_table)

    %% load signal
    fs = gaps_table.fs(i);
    signal = gaps_table.clean{i};

    % write reference wav file for future use
    signal_48 = resample(signal, 48000, fs);
    audiowrite(filename + "_signal.wav", signal_48, 48000);

    for j = 1:length(gaplengths)

        row = row + 1;

        % if the experiment has been stopped, skip what has been already done
        if ~isempty(tables.(methods(end)).PEAQ{row})
            skipped = skipped + 1;
            continue
        end

        % command window output
        signame = string(gaps_table.Properties.RowNames{i});
        str = sprintf("Signal: %s", signame);
        strline = string(repmat('=', 1, strlength(str)));
        fprintf(strline + "\n")
        fprintf(str)
        fprintf("\nGap length: %d ms\n", gaplengths(j))
        fprintf(strline + "\n")

        % estimate the remaining time
        c = datetime("now");
        d = hours(c - cinit); % elapsed time in hours
        fprintf("Elapsed time: %d hours\n", round(d))
        fprintf("Estimated remaining time: %d hours\n\n", ...
            round(d*((rows-skipped)/(row-skipped)-1)))

        %% load mask
        mask = gaps_table.("mask" + num2str(gaplengths(j))){i};
        gapped = signal .* mask;

        %% write meta data to sigtable
        for m = 1:length(methods)
            sigtables.(methods(m)).signal(j) = gaps_table.Properties.RowNames{i};
            sigtables.(methods(m)).fs(j) = gaps_table.fs(i);
            sigtables.(methods(m)).gap(j) = gaplengths(j);
            sigtables.(methods(m)).clean{j} = gaps_table.clean{i};
            sigtables.(methods(m)).mask{j} = gaps_table.("mask" + num2str(gaplengths(j))){i};
            sigtables.(methods(m)).time(j) = 0;
        end

        %% process by segments
        spain = gapped;
        spainmod = gapped;
        
        starts = find(diff(mask) == -1) + 1;
        ends = find(diff(mask) == 1);

        for k = 1:length(starts)

            % A-SPAIN
            [spas, spae] = min_sig_supp_2(...
                w, a, 0, starts(k), ends(k), length(gapped), 1, ...
                offset(starts(k), ends(k), a, 'half'));
            spae = spae - 1;
            SPA.mask = mask(spas:spae);
            SPA.Ls = spae-spas+1;
            t = tic;
            spain(spas:spae) = spain_segmentation(...
                    gapped(spas:spae), ... % degraded signal
                    SPA, ... % model parameters
                    SPA, ... % algorithm parameters
                    signal(spas:spae)); % clean signal (to compute SNR)
            sigtables.(methods(1)).time(j) = sigtables.(methods(1)).time(j) + toc(t);

            % A-SPAIN-MOD
            SPM.mask = mask(spas:spae);
            t = tic;
            restored = a_spain_learned(...
                gapped(spas:spae), ... % degraded signal
                SPM, ... % model parameters
                SPM, ... % algorithm parameters
                signal(spas:spae), ... % clean signal (to compute SNR)
                eye(SPM.M/2+1)); % no learning
            spainmod(spas:spae) = restored(1:(spae-spas+1));
            sigtables.(methods(2)).time(j) = sigtables.(methods(2)).time(j) + toc(t);

        end

        sigtables.(methods(1)).restored{j} = spain;
        sigtables.(methods(2)).restored{j} = spainmod;

        %% compute the metrics
        for m = 1:length(methods)
            iterations = size(sigtables.(methods(m)).restored{j}, 2);
            SDR = NaN(iterations, 1);
            PEMOQ = NaN(iterations, 1);
            PEAQ = NaN(iterations, 1);
            for it = 1:iterations
                solution = sigtables.(methods(m)).restored{j}(:, it);

                if sum(isnan(solution)) > 0
                    break
                end

                % compute SDR
                SDR(it) = snr(signal(~mask), signal(~mask)-solution(~mask));

                % compute PEMO-Q
                [~, ~, PEMOQ(it), ~] = audioqual_silent(signal, solution, fs);
                
                % save the restored signal as wav
                solution_48 = resample(solution, 48000, fs);
                audiowrite(filename + "_solution.wav", solution_48, 48000);

                % compute PEAQ
                PEAQ(it) = PQevalAudio_fn(...
                    filename + "_signal.wav", ...
                    filename + "_solution.wav", 0, length(solution_48));
                delete(filename + "_solution.wav")
            end
            sigtables.(methods(m)).SDR{j} = SDR;
            sigtables.(methods(m)).PEMOQ{j} = PEMOQ;
            sigtables.(methods(m)).PEAQ{j} = PEAQ;

            [M, I] = max(SDR);
            fprintf("%12s, SDR: %.2f dB\n", methods(m), M)
        end

        fprintf("\n")

        %% copy metadata and metrics to tables
        for m = 1:length(methods)
            tables.(methods(m))(row, :) = sigtables.(methods(m))(j, selection);
        end

        S = struct("SPA", SPA, "SPM", SPM, "tables", sigtables);
        save("signals/" + filename + "_" + signame + ".mat", "-struct", "S", "-v7.3")
        save(filename + ".mat", "SPA", "SPM", "tables", "-v7.3")
        
    end

    delete(filename + "_signal.wav")

end
