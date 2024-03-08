function y = arinpaint(x, maxlen, order, method)
% ARINPAINT Inpainting based on a forward and backward linear prediction

% find the missing samples, identified as NaN
mask = ~isnan(x);
s = find(~mask, 1, "first");
f = find(~mask, 1, "last");
h = f - s + 1;

% manage inputs
if nargin < 2
    presig = x(1:s-1);
    postsig = x(f+1:end);
else
    presig = x(max(1, s-maxlen):s-1);
    postsig = x(f+1:min(length(x), f+maxlen));
end
if nargin < 3
    order = max(length(presig), length(postsig)) - 1;
end

% forward prediction
premean = mean(presig);
presig = presig - premean;
if strcmpi(method, "lpc")
    af = lpc(presig, order);
else
    af = arburg(presig, order);
end
zf = filtic(1, af, flip(presig(end-order+1:end)));
prediction = filter(1, af, zeros(1, h), zf)';

% backward prediction
postsig = flip(postsig);
postmean = mean(postsig);
postsig = postsig - postmean;
if strcmpi(method, "lpc")
    ab = lpc(postsig, order);
else
    ab = arburg(postsig, order);
end
zb = filtic(1, ab, flip(postsig(end-order+1:end)));
postdiction = flip(filter(1, ab, zeros(1, h), zb)');

% composition
y = x;
t = linspace(0, pi/2, h)';
wts = cos(t).^2;
% wts = (h:-1:1)'/(h+1); % this wights are used in fillgaps
y(s:f) = wts .* (prediction + premean) + (1 - wts) .* (postdiction + postmean);

end

