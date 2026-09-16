function xrec = Etter(x, p, w, optional)

    arguments
        x (:, 1) {mustBeNumeric}
        p (1, :) {mustBeNumeric}
        w (1, 1) {mustBeNumeric}
        optional.l (1, :) {mustBeInteger, mustBeNonnegative} = [] % gap starts in samples
        optional.M (1, :) {mustBeInteger, mustBeNonnegative} = [] % gap lengths in samples
        optional.estimateFunc (1, 1) string {mustBeMember(...
            optional.estimateFunc, ["lpc", "arburg"])} = "lpc"
        optional.mode (1, 1) string {mustBeMember(...
            optional.mode, ["inpainting", "plc"])} = "inpainting"
    end
    
    if isempty(optional.l) || isempty(optional.M)
        [optional.M, optional.l] = FindGaps(x);
    end
    
    for i = 1:length(optional.M)
        xrec = CalculateEtter(x, p(i), optional.l(i), ...
            optional.M(i), w, optional.estimateFunc, optional.mode);
    end
end

function xrec = CalculateEtter(x, p, l, M, w, estimateFunc, mode)
    xrec = x;

    %% Define A and B
    wl = max(1, l - w);
    wr = min(length(x), l + M + w - 1);

    switch estimateFunc
        case "lpc"
            a = lpc(x(wl:l-1), p);
            if mode == "inpainting"
                b = lpc(fliplr(x(l+M:wr)), p);
            end
        case "arburg"
            a = arburg(x(wl:l-1), p);
            if mode == "inpainting"
                b = arburg(fliplr(x(l+M:wr)), p);
            end
    end

    colA = zeros(M, 1);
    lenA = min(M, p + 1);
    colA(1:lenA) = a(1:lenA);
    A = toeplitz(colA, [a(1), zeros(1, M - 1)]);

    if mode == "inpainting"
        rowB = zeros(1, M);
        lenB = min(M, p + 1);
        rowB(1:lenB) = b(1:lenB);
        B = toeplitz([b(1); zeros(M - 1, 1)], rowB);
    end

    %% Define L and R
    L = zeros(M, p + 1);
    for i = 1:M
        for j = (i + 1):(p + 1)
            L(i, j) = x(l + i - j);
        end
    end
    if mode == "inpainting"
        R = zeros(M, p + 1);
        for i = 1:M
            for j = 1:(p + 1)
                if j > M - i + 1
                    R(i, j) = x(l + i + j - 2);
                end
            end
        end
    end

    %% Compute matrix D and vector y
    if mode == "inpainting"
        % D = A'*A + B'*B;
        % y = -A'*L*a'-B'*R*b';

        % [OM] brackets can help
        % y = -A'*(L*a')-B'*(R*b');

        % [OM] Fourier can help
        N = M+p;
        AA = ifft(fft(flip(a(:)), N) .* fft(A, N));
        AA = AA(p+1:end, :);

        BB = ifft(fft(b(:), N) .* fft(B, N));
        BB = BB(1:end-p, :);
        D = AA + BB;

        y1 = -ifft(fft(flip(a(:)), N) .* fft(L*a', N));
        y1 = y1(p+1:end);
        y2 = -ifft(fft(b(:), N) .* fft(R*b', N));
        y2 = y2(1:end-p);
        y = y1 + y2;

        % check
        % norm(A'*A-AA)
        % norm(B'*B-BB)
        % norm(y1-(-A'*L*a'))
        % norm(y2-(-B'*R*b'))
    end

    %% Solve the equation and compute the missing signal
    if mode == "inpainting"
        xn = D \ y;
    elseif mode == "plc"
        xn = A \ (-L * a(:));
    end
    xrec(l:l + M - 1) = xn.';

end

function [M, l] = FindGaps(x)
    zeroMask = isnan(x);
    edges = diff([0; zeroMask; 0]);
    l_tmp = find(edges == 1); % gap starts
    gaps_end = find(edges == -1) - 1; % gap ends
    M_tmp = gaps_end - l_tmp + 1; % gap lengths

    % Remove overly short gaps/zero values
    valid_idx = M_tmp > 2;
    l = l_tmp(valid_idx);
    M = M_tmp(valid_idx);
end