%% Function: Simulate signal detection
%

function [Nhit, Nfa, r]=SimEVSD(param, Ntarget, Nlure, Nsubj, Nc)
% param: [d-prime, c]
dp=param(1); % d-prime (sensitivity)
c=param(2:Nc+1); % criterion
p_hit=1-normcdf(c, dp, 1); % hit rate
p_fa=1-normcdf(c, 0, 1); % false-alarm rate
Nhit=zeros(Nsubj,Nc);
Nfa=Nhit;
r=rng; % record random seeds
for subj=1:Nsubj
    for c_id=1:Nc
        Nhit(subj, c_id)=binornd(Ntarget(c_id), p_hit(c_id)); % hit trials
        Nfa(subj, c_id)=binornd(Nlure(c_id), p_fa(c_id)); % false-alarm trials
    end
end
end