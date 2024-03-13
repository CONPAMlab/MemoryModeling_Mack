%% Get drift functions
%
%

function df = GetDF(res)

w = res.wts;

n          = 100;
kBase      = sd2k((2*pi/numel(w)));
dx = 2*pi/n; %bin width
xe = linspace(-pi,pi,n+1); %bin edges
xc = (xe(1:n) + dx/2)'; %bin centers

muVM = -pi:(2 * pi / numel(w)):(pi - 2 * pi / numel(w));
df   = sum(w' .* (kBase * sin(muVM - xc) .* exp(kBase * cos(muVM - xc))),2);
% df   = sum(w' .* (kBase * sin(muVM - xc) .* exp(kBase * cos(muVM - xc))) ...
%     ./sum(2 * pi * besseli(0,kBase)),2);
df   = df / max(abs(df));

end

function K = sd2k (S)
% SD2K (S)
%   Returns the Von Mises concentration parameter corresponding to
%   standard deviation S of a wrapped normal distribution.
%
%   Ref: Topics in Circular Statistics, S. R. Jammalamadaka & A. Sengupta
%
%   --> www.paulbays.com

R = exp(-S.^2/2);
K = 1./(R.^3 - 4 * R.^2 + 3 * R);
K(R < 0.85) = -0.4 + 1.39 * R(R < 0.85) + 0.43./(1 - R(R < 0.85));
K(R < 0.53) = 2 * R(R < 0.53) + R(R < 0.53).^3 + (5 * R(R < 0.53).^5)/6;
end
