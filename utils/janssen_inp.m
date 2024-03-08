function [restored, objective, times] = janssen_inp(signal, p, maxit, varargin)
% janssen performs inpainting using the Janssen algorithm [1]
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
%   signal        the input (degraded) signal, missing samples identified
%                 with NaN values
%   p             order of the AR model
%   maxit         number of iterations of the whole Janssen algorithm
%   varargin      name-value pairs
%                 "method" ("lpc")  switch between AR model estimators,
%                                   accepted are options "lpc" and "arburg"
%                 "saveall" (false) save the solution during iterations
%                 "verbose" (false) print current iteration
%
% output arguments
%   restored      the solution; if saveall is true, restored is of size
%                 length(signal) x maxit, or length(signal) x 1 otherwise
%   objective     values of the objective function during iterations,
%                 computed using the forward error only, i.e., it only
%                 corresponds to the true objective with method set to
%                 "lpc"
%   times         cumulative computation time during iterations
%
% Date: 23/02/2024
% By Ondrej Mokry
% Brno University of Technology
% Contact: ondrej.mokry@vut.cz

%% parse the inputs
% create the parser
pars = inputParser;
pars.KeepUnmatched = true;

% add optional name-value pairs
addParameter(pars, "method", "lpc")
addParameter(pars, "saveall", false)
addParameter(pars, "verbose", false)

% parse
parse(pars, varargin{:})

% save the parsed results to nice variables
method  = pars.Results.method;
saveall = pars.Results.saveall;
verbose = pars.Results.verbose;

%% initialization
mask = ~isnan(signal);
solution = signal;
solution(~mask) = 0;
N = length(signal);
if saveall
    restored = NaN(N, maxit);
end
if nargout > 1
    Q = @(x, c) 0.5*norm(fft(c, N+p).*fft(x, N+p))^2 / (N+p);
    objective = NaN(maxit, 1);
end

% prepare some matrices
indmiss = find(~mask);
indobs  = find(mask);
IAA     = abs(repmat(indmiss, 1, N)-repmat(1:N, length(indmiss), 1));
IAA1    = IAA <= p;

% if desired, start the timer
if nargout > 2
    times = NaN(maxit, 1);
    tic
end

%% main iteration
if verbose
    str = "";
end
for i = 1:maxit

    if verbose
        fprintf(repmat('\b', 1, strlength(str)))
        str = sprintf("iteration %d of %d", i, maxit);
        fprintf(str)
    end

    % AR model estimation
    if strcmpi(method, "lpc")
        coef = lpc(solution, p)'; 
    else
        coef = arburg(solution, p)';
    end

    % signal estimation
    AA       = zeros(size(IAA));
    b        = coef'*hankel(coef', [coef(end), zeros(1, p)]);
    AA(IAA1) = b(IAA(IAA1)+1);
    [R, er]  = chol(AA(:, indmiss));
    if er
        break
    else
        solution(~mask) = -R\(R'\(AA(:, indobs)*signal(indobs)));
    end

    % update the solution
    if saveall
        restored(:, i) = solution;
    else
        restored = solution;
    end
    
	% compute the objective value
    if nargout > 1
        objective(i) = Q(solution, coef);
    end

    % update the elapsed time
    if nargout > 2
        times(i) = toc;
    end

end

if verbose
    fprintf(repmat('\b', 1, strlength(str)))
end

end