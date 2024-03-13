%% Function: Simulate double-threshold dual-process signal detection
%

function [data, r]=SimDPSD(param, Ntarget, Nlure, Nsubj, Nc)
% param: [d-prime, c, sigma, Ro, Rn]
dp=param(1); % d-prime (sensitivity)
c=param(2:Nc+1); % criterion
sigma=param(Nc+2); % SD for target density
Ro=param(Nc+3); % p(always report target)
Rn=param(Nc+4); % p(always report noise)
p_hit=Ro+(1-Ro)*(1-normcdf(c, dp, sigma)); % hit rate
p_fa=(1-Rn)*(1-normcdf(c, 0, 1)); % false-alarm rate
r=rng; % record random seeds
data = cell(Nsubj, 1);
for subj=1:Nsubj
    Nhit=zeros(1,Nc);
    Nfa=Nhit;
    for c_id=1:Nc
        Nhit(c_id)=binornd(Ntarget, p_hit(c_id)); % hit trials
        Nfa(c_id)=binornd(Nlure, p_fa(c_id)); % false-alarm trials
    end
    data{subj}.Nhit = [Ntarget Nhit];
    data{subj}.Nfa = [Nlure Nfa];
end
end
% for i=1:100
%
%     a(i)=binornd(10000,.5)
% end