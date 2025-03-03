clear
clc
close all

gapnum = 10;
gaplengths = 10:10:80;

%% load data
fold = "IRMAS";
files = struct2table(dir(fold));
files = files(~files.isdir, :);
data.fs = 44100;
names = strings(height(files), 1);
for i = 1:height(files)
    names(i) = sprintf("irmas%02d", i);

    % load signal
    signal = audioread(fullfile(fold, files.name{i}));
    
    % convert to mono
    signal = mean(signal, 2);

    % take 7 seconds around the middle of the signal
    indices = round(length(signal)/2 - 3.5*data.fs) + (1:7*data.fs);
    signal = fade(signal(indices), round(0.1*data.fs));

    % write to the struct data
    data.(names(i)) = signal;
end
clear indices signal

%% prepare table
tbl = table('Size', [length(names), 11], ...
    'VariableTypes', ...
    ["string", "double", "cell", "cell", "cell", "cell", "cell", ...
    "cell", "cell", "cell", "cell"], ...
    'VariableNames', ...
    ["filename", "fs", "clean", "mask10", "mask20", "mask30", "mask40", ...
    "mask50", "mask60", "mask70", "mask80"], ...
    'RowNames', names);

signal = [];
gap = [];
shift = [];
shift = round(shift * data.fs);
shifts = table(signal', gap', shift', 'VariableNames', ["signal", "gap", "shift"]);
clear signal gap shift

%% fill table
for i = 1:length(names)
    tbl.fs(i) = data.fs;
    tbl.clean{i} = data.(names{i});
    tbl.filename(i) = string(files.name{i});
    for gaplength = gaplengths
        % generate mask
        [~, mask] = makegaps(tbl.clean{i}, tbl.fs(i), gapnum, gaplength, "w", 8192);
        starts = find(diff(mask) < 0);
        ends = find(diff(mask) > 0);

        % modify
        rows = find(shifts.signal == i);
        for r = rows
            starts(shifts.gap(r)) = starts(shifts.gap(r)) + shifts.shift(r);
            ends(shifts.gap(r)) = ends(shifts.gap(r)) + shifts.shift(r);
        end
        diffm = zeros(length(mask)-1, 1);
        diffm(starts) = -1;
        diffm(ends) = 1;
        mask = cumsum([1; diffm]);
        mask = logical(mask);

        % find minimal distance between gaps
        dists = zeros(gapnum+1, 1);
        dists(1) = starts(1)-1; % number of clean samples from start
        dists(end) = length(tbl.clean{i})-ends(end); % number of clean samples from end
        for j = 1:gapnum-1
            dists(1+j) = starts(j+1)-ends(j)-1;
        end
        fprintf("%s, gap length %d ms: min gap distance %.1f ms\n", ...
            names{i}, gaplength, min(dists)/tbl.fs(i)*1000)

        % write to table
        tbl.("mask" + gaplength){i} = mask;
    end
end

%% save as table
gaps_table = tbl;
save("gaps_table.mat", "gaps_table")
clear gaps_table

% read example for signal number 2, gap length 20 ms:
%
% load("gaps_table.mat")
% clean_signal = gaps_table.("clean"){"irmas02"};
% mask = gaps_table.("mask20"){"irmas02"};
%
% 1st argument: "mask" + length in miliseconds
% 2nd argument: signal name

%% save as struct
gaps_struct = table2struct(tbl);
save("gaps_struct.mat", "gaps_struct")
clear gaps_struct

% read example for violin signal, gap length 20 ms:
%
% load("gaps_struct.mat")
% clean_signal = gaps_struct(1).("clean");
% mask = gaps_struct(1).("mask20");
%
% 1st argument: signal number (starting from 1)
% 2nd argument: "mask" + length in miliseconds

% read example in python for violin signal, gap length 20 ms:
%
% import scipy.io as sio
% gaps_sio = sio.loadmat("gaps_struct.mat")
% clean_signal = gaps_sio["gaps_struct"]["clean"][0, 0]
% mask = gaps_sio["gaps_struct"]["mask20"][0, 0]
%
% 1st argument: "gaps_struct"
% 2nd argument: "mask" + length in miliseconds
% 3rd argument: signal id (starting from 0)
% 4th argument: 0

%% plot
% for i = 1:height(tbl)
% 
%     figure("Name", names{i})
%     tls = tiledlayout("flow");
%     title(tls, names{i}, "Interpreter", "none")
%     time = (0:length(tbl.clean{i})-1)/tbl.fs(i);
% 
%     for gaplength = gaplengths
% 
%         nexttile
%         plot(time, tbl.clean{i}, time, max(tbl.clean{i}) * tbl.("mask" + gaplength){i})
%         axis tight
%         title("gap " + num2str(gaplength) + " ms")
%     end
% 
% end

function y = fade(x, len)
    
    y = x;
    cosinus = cos(linspace(0, pi/2, len)').^2;
    y(1:len, :) = y(1:len, :) .* (1-cosinus);
    y(end-len+1:end, :) = y(end-len+1:end, :) .* cosinus;

end