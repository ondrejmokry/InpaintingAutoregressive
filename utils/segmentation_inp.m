function restored = segmentation_inp(signal, p, maxit, varargin)
% segmentation performs segment-wise audio inpainting using the Janssen
% algorithm [1]
%
% implementation of the signal estimation step is taken from the Audio
% Inpainting Toolbox as described in [2]
%
% [1] A. Janssen, R. Veldhuis and L. Vries, "Adaptive interpolation of
%     discrete-time signals that can be modeled as autoregressive
%     processes, " in IEEE Transactions on Acoustics, Speech, and Signal
%     Processing, vol. 34, no. 2, pp. 317-330, 1986, doi:
%     10.1109/TASSP.1986.1164824.
% [2] A. Adler, V. Emiya, M. G. Jafari, M. Elad, R. Gribonval and M. D.
%     Plumbley, "Audio Inpainting, " in IEEE Transactions on Audio, Speech, 
%     and Language Processing, vol. 20, no. 3, pp. 922-932, 2012, doi:
%     10.1109/TASL.2011.2168211.
%
% input arguments
%   signal        the input (degraded) signal
%   p             order of the AR model
%   maxit         number of iterations of the whole Janssen algorithm
%   varargin      name-value pairs
%                 'wtype' ('hann')    window shape
%                 'w' (4096)          window length
%                 'a' (2048)          window shift
%                 'verbose' (true)    print number of the segment being
%                                     processed
%
% output arguments
%   restored      the solution; if saveall is true, restored is of size
%                 length(signal) x maxit, or length(signal) x 1 otherwise
%
% Date: 23/02/2024
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

%% parse the inputs for overlap-add
% create the parser
pars = inputParser;
pars.KeepUnmatched = true;

% add optional name-value pairs
addParameter(pars, 'wtype', 'hann')
addParameter(pars, 'w', 4096)
addParameter(pars, 'a', 1024)
addParameter(pars, 'saveall', false)
addParameter(pars, 'verbose', true)

% parse
parse(pars, varargin{:})

% save the parsed results to nice variables
wtype   = pars.Results.wtype;
wtype   = char(wtype);
w       = pars.Results.w;
a       = pars.Results.a;
saveall = pars.Results.saveall;
verbose = pars.Results.verbose;

%% padding the signal
L    = ceil(length(signal)/a)*a + (ceil(w/a)-1)*a;
S    = L/a; % number of signal segments
data = [signal; zeros(L-length(signal), 1)];

%% initializing the solution array
if saveall
    mrestored = zeros(w, maxit, S);
else
    mrestored = zeros(w, S);
end

%% construction of analysis and synthesis windows
if strcmpi(wtype, 'rect')
    gana = ones(w, 1);
    gsyn = gabwin('hann', a, w, L);
    gsyn = fftshift(gsyn);
    gsyn = normalize(gsyn, 'peak');
elseif strcmpi(wtype, 'tukey')
    gana = tukeywin(w);
    gsyn = gana; % this will be compensated later
else
    g    = gabwin(wtype, a, w, L);
    gana = normalize(g, 'peak'); % peak-normalization of the analysis window
    gana = fftshift(gana);
    gsyn = gabdual(gana, a, w)*w; % computing the synthesis window
end

%% segment settings
mdata = NaN(w, S);
for s = 1:S
    % defining the indices of the current block
    indices = 1 + (s-1)*a - floor(w/2) : (s-1)*a + ceil(w/2);
    indices = 1 + mod(indices-1, L);
    
    % defining the segment data and masks
    mdata(:, s) = data(indices) .* gana;
end

%% segment processing via parfor
if saveall
    parfor s = 1:S
        if verbose
            fprintf('Processing segment %d of %d...\n', s, S)
        end
        if sum(isnan(mdata(:, s))) == w % all samples are missing
            mrestored(:, :, s) = 0;
        elseif sum(isnan(mdata(:, s))) == 0 % no samples are missing
            mrestored(:, :, s) = repmat(mdata(:, s), [1, maxit]);
        else
            mrestored(:, :, s) = janssen_inp(mdata(:, s), p, maxit, varargin{:}); %#ok<PFBNS>
        end
    end
else
    parfor s = 1:S
        if verbose
            fprintf('Processing segment %d of %d...\n', s, S)
        end
        if sum(isnan(mdata(:, s))) == w % all samples are missing
            mrestored(:, s) = 0;
        elseif sum(isnan(mdata(:, s))) == 0 % no samples are missing
            mrestored(:, s) = mdata(:, s);
        else
            mrestored(:, s) = janssen_inp(mdata(:, s), p, maxit, varargin{:}); %#ok<PFBNS>
        end
    end
end

%% overlap-add
rescale = zeros(L, 1);
if saveall
    restored = zeros(L, maxit);
else
    restored = zeros(L, 1);
end
for s = 1:S
    indices = 1 + (s-1)*a - floor(w/2) : (s-1)*a + ceil(w/2);
    indices = 1 + mod(indices-1, L);
    rescale(indices) = rescale(indices) + gana.*gsyn;
    if saveall
        restored(indices, :) = restored(indices, :) + mrestored(:, :, s).*repmat(gsyn, 1, maxit);
    else
        restored(indices) = restored(indices) + mrestored(:, s).*gsyn;
    end
end

%% rescale
restored = restored ./ rescale;

%% cropping the solution to the original length
restored = restored(1:length(signal), :);

end