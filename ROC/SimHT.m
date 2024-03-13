%% Function: Simulate HT Models
%

function [data, r]=SimHT(param, Ntarget, Nlure, Nsubj, Ng)
% param: [c, Ro, Rn]
g=param(1:Ng); % criterion
Ro=param(Ng+1); % p(always report target)
Rn=param(Ng+2); % p(always report noise)
p_hit=Ro+(1-Ro)*g; % hit rate
p_fa=(1-Rn)*g; % false-alarm rate
r=rng; % record random seeds
data = cell(Nsubj, 1);
for subj=1:Nsubj
    Nhit=zeros(1,Ng);
    Nfa=Nhit;
    for c_id=1:Ng
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