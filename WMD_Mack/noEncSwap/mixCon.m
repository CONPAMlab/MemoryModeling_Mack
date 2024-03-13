function [c, ceq] = mixCon(x,D,SS)
%mixing weights must be between 0 and 1

nSS = numel(unique(SS));
puB    = x((2 * nSS+1):(2 * nSS+4));

pU  = puB(1) + SS .* puB(2) + D .* puB(3) + SS .* D .* puB(4);
% pS  = psB(1) + SS .* psB(2) + D .* psB(3) + SS .* D .* psB(4);
% pS(SS==1) = 0;
pM = (1-pU);

c(1) = max(pM) - 1; %pM <= 1 
c(2) = -min(pM); %pM >= 0
c(3) = max(pU) - 1; 
c(4) = -min(pU); 

ceq = [];