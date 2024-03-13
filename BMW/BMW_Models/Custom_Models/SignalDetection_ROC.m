%% Function: Fit double-threshold dual-process signal detection
%

function Output=SignalDetection_ROC(param, Data, Input)
% Specify parameters
dp=param(1); % d-prime (sensitivity)
Nc=length(Data.Ntarget); % number of bins-1
c=param(2:1+Nc); % response criterion
Nparam=Nc+1;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
if any(strcmp(Input.Variants,'Unequal Variance'))
    Nparam=Nparam+1;
    sigma=param(Nparam); % SD for target density
else
    sigma=1;
end
if any(strcmp(Input.Variants,'Recollection Threshold'))
    Nparam=Nparam+1;
    Ro=param(Nparam); % p(always report target)
else
    Ro=0;
end
if any(strcmp(Input.Variants,'Newness Threshold'))
    Nparam=Nparam+1;
    Rn=param(Nparam); % p(always report target)
else
    Rn=0;
end
if any(strcmp(Input.Variants,'Recognition Rate'))
    Nparam=Nparam+1;
    Pm=param(Nparam); % p(always report target)
else
    Pm=1;
end
Ntarget=Data.Nhit(1);
Nlure=Data.Nfa(1);
Nhit=Data.Nhit(2:end);
Nfa=Data.Nfa(2:end);
if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end

LLH=0;
p_hit=zeros(1,Nc);
p_fa=zeros(1,Nc);
for i=1:Nc
    p_hit(i)=Pm*(Ro+(1-Ro)*(1-normcdf(c(i), dp, sigma)))+(1-Pm)*(1-normcdf(c(i), 0, 1)); % hit rate
    p_fa(i)=(1-Rn)*(1-normcdf(c(i), 0, 1)); % false-alarm rate
    LLH_target=-log(binopdf(Nhit(i), Ntarget, p_hit(i)));
    LLH_lure=-log(binopdf(Nfa(i), Nlure, p_fa(i)));
    LLH=LLH+LLH_target+LLH_lure; % - log likelihood
end

% Posterior
LP=-log(Prior)+LLH; % likelihood*prior

if LP==Inf || isnan(LP)
    LP=realmax('double'); % Output should be a real value
    LLH=LP;
end

% PPD
p_LH=zeros(sum(Ntarget)+sum(Nlure),1);
N_LH=0;

for i=1:Nc
    p_LH(N_LH+1:N_LH+Ntarget(i)+Nlure(i),1)=[p_hit(i)*ones(Nhit(i),1); (1-p_hit(i))*ones(Ntarget(i)-Nhit(i),1); p_fa(i)*ones(Nfa(i),1); (1-p_fa(i))*ones(Nlure(i)-Nfa(i),1)];
    N_LH=N_LH+Ntarget(i)+Nlure(i);
end

% Decide output
if strcmp(Input.Output,'LP')
    Output=LP;
elseif strcmp(Input.Output,'LLH')
    Output=LLH;
elseif strcmp(Input.Output,'Prior')
    Output=Prior;
elseif strcmp(Input.Output,'LPPD')
    Output=log(p_LH);
elseif strcmp(Input.Output,'All')
    Output.LP=LP;
    Output.LLH=LLH;
    Output.Prior=Prior;
    Output.LPPD=log(p_LH);
end

end

% Define prior
function p=prior(param, Input)

p=1;

end
