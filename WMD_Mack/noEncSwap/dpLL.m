function ll = dpLL(R,T,D,SS,x)
           
nSS = numel(unique(SS));
sigmaM = x((0 * nSS + 1):(1 * nSS));
betaM  = x((1 * nSS + 1):(2 * nSS));
puB    = x((2 * nSS + 1):(2 * nSS + 4));
w      = x((2 * nSS + 5):(end-1));
sigmaR = x(end);

% [P, xc] = dp(T, NT, D, SS, sigmaM, sigmaE, betaM, betaE, puB, psB, w, sigmaR);
[P, xc] = dp(T, D, SS, sigmaM, betaM, puB, w, sigmaR);
[~, i] = min(abs(xc - R'),[],1); 
ind = sub2ind(size(P),i,1:numel(T)); 
p = P(ind);
ll = sum(-log(p));
if ~isreal(ll); ll = inf; end

