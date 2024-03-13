%% Play around with different basis for drifting function
%
%
close all

w = [20 10 20 10 1 1 50];
w = (w./sum(w))';
n = 100;
kBase = sd2k((2*pi/numel(w)));

dx = 2*pi/n; %bin width
xe = linspace(-pi,pi,n+1); %bin edges
xc = (xe(1:n) + dx/2)'; %bin centers

%drift function
muVM = -pi:(2 * pi / numel(w)):(pi - 2 * pi / numel(w)); % w = weights of the basis

muVM = [11 42 73 104 130 148 170].*2/180*pi-pi;

% df   = sum(w' .* (kBase * sin(muVM - xc) .* exp(kBase * cos(muVM - xc))) ...
%     ./sum(2 * pi * besseli(0,kBase)),2);
df   = sum(w' .* (kBase * sin(muVM - xc) .* exp(kBase * cos(muVM - xc))),2);
% df   = sum(w' .* (exp(kBase * cos(muVM - xc)))...
%         ./sum(2 * pi * besseli(0,kBase)),2);
% df   = sum(w' .* (exp(kBase * cos(muVM - xc))),2);
df = sum(w' .* (2./(1+exp(kBase*(muVM-xc)))-1), 2);

df   = df / max(abs(df));

plot(xc, (kBase * sin(muVM - xc) .* exp(kBase * cos(muVM - xc)))./max(kBase * sin(muVM - xc) .* exp(kBase * cos(muVM - xc))),'k:')
plot(xc, 2./(1+exp(kBase*(muVM-xc)))-1, 'k:')
hold on
plot(xc, df)
plot(muVM, w, 'o')
ylim([-1,1])
xlim([-pi,pi])

% haha = -pi:0.01:pi;
% plot(kBase*sin(haha).*exp(kBase*cos(haha)))