clear
clc
close all

paperDir = fileparts(mfilename("fullpath"));
repoDir = fileparts(fileparts(paperDir));
resultsDir = fullfile(paperDir, "results 2");
signalsDir = fullfile(paperDir, "signals");
if ~isfolder(resultsDir)
    mkdir(resultsDir)
end
if ~isfolder(signalsDir)
    mkdir(signalsDir)
end
addpath(genpath(fullfile(repoDir, "utils")))

fileid = 0;
for p = [256, 512, 1024, 2048, 3072, 4096]

    for method = ["arburg", "lpc"]

        fileid = fileid + 1;

        % future file name
        filename = fullfile(resultsDir, sprintf("results_%02d", fileid));

        %% load signals and masks
        if contains(resultsDir, "irmas")
            load(fullfile(repoDir, "irmas_gaps_table.mat"))
        else
            load(fullfile(repoDir, "gaps_table.mat"))
        end
        gaplengths = 10:10:80; % this corresponds to what is saved in gaps_table
        
        %% set params
        w = 4096;
        a = w/4;
        maxit = 20;

        %% skip invalid combination
        if w == 4096 && p == 4096
            continue
        end
        
        % if the experiment has been stopped, load the data
        if isfile(filename + ".mat")
            s = load(filename + ".mat");

            % check compatibility; load table only if it is compatible
            if s.a ~= a || s.method ~= method || s.p ~= p || s.w ~= w
                warning("The saved table does not fit the script, sorry.")
            else
                tables = s.tables;
            end
            clear s
        else
            tables = struct;
        end
        skipped = 0;

        %% prepare table for saving everything
        methods = ["etter", "etter_plc"];
        types = ["string", "double", "double", "cell", "cell", ...
            "string", "double", "double", "double", "double", ...
            "cell", "double", "cell", "cell", "cell"];
        names = ["signal", "fs", "gap", "clean", "mask", ...
            "method", "p", "w", "a", "maxit", ...
            "restored", "time", "SDR", "PEMOQ", "PEAQ"];
        units = ["", "Hz", "ms", "", "", "", "", "", "", "", "", "s", "dB", "", ""];
        selection = [1, 2, 3, 6:10, 12:15]; % all except clean, mask and restored signals
        rows = height(gaps_table) * length(gaplengths);
        for m = 1:length(methods)
            sigtables.(methods(m)) = table('Size', [length(gaplengths), length(types)], ...
                'VariableTypes', types, 'VariableNames', names);
            sigtables.(methods(m)).Properties.VariableUnits = units;
            if ~isfield(tables, methods(m))
                tables.(methods(m)) = table('Size', [rows, length(selection)], ...
                    'VariableTypes', types(selection), 'VariableNames', names(selection));
                tables.(methods(m)).Properties.VariableUnits = units(selection);
            end
        end

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
                str_1 = sprintf("Signal: %s", signame);
                str_2 = sprintf("Gap length: %d ms", gaplengths(j));
                strline = string(repmat('=', 1, max(strlength(str_1), strlength(str_2))));
                fprintf(strline + "\n")
                fprintf(str_1 + "\n" + str_2 + "\n")
                fprintf(strline + "\n")
        
                % estimate the remaining time
                c = datetime("now");
                d = hours(c - cinit); % elapsed time in hours
                fprintf("Elapsed time: %d hours\n", round(d))
                fprintf("Estimated remaining time: %d hours\n\n", ...
                    round(d*((rows-skipped)/(row-skipped)-1)))
        
                %% load mask
                mask = gaps_table.("mask" + num2str(gaplengths(j))){i};
                gapped = signal;
                gapped(~mask) = NaN;
        
                %% write meta data to sigtable
                for m = 1:length(methods)
                    sigtables.(methods(m)).signal(j) = gaps_table.Properties.RowNames{i};
                    sigtables.(methods(m)).fs(j) = gaps_table.fs(i);
                    sigtables.(methods(m)).gap(j) = gaplengths(j);
                    sigtables.(methods(m)).clean{j} = gaps_table.clean{i};
                    sigtables.(methods(m)).mask{j} = gaps_table.("mask" + num2str(gaplengths(j))){i};
                    sigtables.(methods(m)).method(j) = method;
                    sigtables.(methods(m)).p(j) = p;
                    sigtables.(methods(m)).w(j) = w;
                    sigtables.(methods(m)).a(j) = a;
                    sigtables.(methods(m)).maxit(j) = maxit;
                    sigtables.(methods(m)).time(j) = 0;
                end
        
                %% process by segments
                etter = gapped;
                etter_plc = gapped;

                starts = find(diff(mask) == -1) + 1;
                ends = find(diff(mask) == 1);
        
                for k = 1:length(starts)

                    segs = starts(k) - w;
                    sege = ends(k) + w;

                    % Etter
                    t = tic;
                    etter(segs:sege) = Etter(gapped(segs:sege), p, w, estimateFunc=method, mode="inpainting");
                    sigtables.(methods(1)).time(j) = sigtables.(methods(1)).time(j) + toc(t);

                    % Etter PLC
                    t = tic;
                    etter_plc(segs:sege) = Etter(gapped(segs:sege), p, w, estimateFunc=method, mode="plc");
                    sigtables.(methods(2)).time(j) = sigtables.(methods(2)).time(j) + toc(t);
                end
        
                sigtables.(methods(1)).restored{j} = etter;
                sigtables.(methods(2)).restored{j} = etter_plc;
        
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
                    if isscalar(SDR)
                        fprintf("%13s, SDR: %.2f dB\n", methods(m), M)
                    else
                        fprintf("%13s, peak SDR: %.2f dB (iteration %d of %d), end SDR: %.2f dB\n", ...
                            methods(m), M, I, maxit, SDR(end))
                    end
                end
        
                fprintf("\n")
        
                %% copy metadata and metrics to tables
                for m = 1:length(methods)
                    tables.(methods(m))(row, :) = sigtables.(methods(m))(j, selection);
                end
        
                S = struct("a", a, "maxit", maxit, "method", method, "p", p, "tables", sigtables, "w", w);
                % save("signals/" + filename + "_" + signame + ".mat", "-struct", "S", "-v7.3")
                save(filename + ".mat", "a", "maxit", "method", "p", "tables", "w", "-v7.3")
                
            end
        
            delete(filename + "_signal.wav")
        
        end
    end
end