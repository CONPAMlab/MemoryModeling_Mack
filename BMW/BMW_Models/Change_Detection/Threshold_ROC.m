%% Function: Fit Threshold ROC models
%

function Output=Threshold_ROC(param, Data, Input)
% Specify parameters
Nc=length(Data.Nhit)-1; % number of bins-1
c=param(1:Nc); % response criterion
Ro=param(Nc+1); % Oldness threshold
Rn = 0; % Newness threshold
Nparam=Nc+1;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
if any(strcmp(Input.Variants,'Double Equal Threshold'))
    Rn=Ro;
end
if any(strcmp(Input.Variants,'Newness Threshold'))
    Nparam=Nparam+1;
    Rn=param(Nparam); % p(always report target)
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
    p_hit(i)=Ro+(1-Ro)*c(i); % hit rate
    p_fa(i)=(1-Rn)*c(i); % false-alarm rate
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
    p_LH(N_LH+1:N_LH+Ntarget+Nlure,1)=[p_hit(i)*ones(Nhit(i),1); (1-p_hit(i))*ones(Ntarget-Nhit(i),1); p_fa(i)*ones(Nfa(i),1); (1-p_fa(i))*ones(Nlure-Nfa(i),1)];
    N_LH=N_LH+Ntarget+Nlure;
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
% 
% % Specify parameters
% K=param(1); % Capacity
% % weibull prior for capacity
% p0(1)=wblpdf(K,3.5,3); % given that K is ofter 3~4
% kappa_1=param(2); % Unit resource
% % Gamma prior for unit resource
% p0(2)=gampdf(kappa_1,3,5);
% Nparam=2;
% 
% if ~isfield(Input,'Variants') % No Variants
%     Input.Variants={};
% end
% if any(strcmp(Input.Variants,'ResponseNoise'))
%     kappa_r=param(3); % Response variability
%     % Gamma prior for response noise
%     p0(3)=gampdf(kappa_r,3,5);
% end
% if any(strcmp(Input.Variants,'Bias'))
%     Nparam=Nparam+1;
%     bias=param(Nparam); % Mean bias
%     % Gaussian prior for bias
%     p0(Nparam)=normpdf(bias, 0, 1);
% end
% if any(strcmp(Input.Variants,'BiasF'))
%     Nparam=Nparam+1;
%     biasF=param(Nparam); % Fluctuation of bias
%     % Gaussian prior for the fluctuation of bias
%     p0(Nparam)=normpdf(biasF, 0, 5);
% end
% if any(strcmp(Input.Variants,'PrecF'))
%     Nparam=Nparam+1;
%     precF=param(Nparam); % Fluctuation of precision
%     % Gaussian prior for the fluctuation of precision
%     p0(Nparam)=normpdf(precF, 0, 1);
% end
% if any(strcmp(Input.Variants,'Swap'))
%     Nparam=Nparam+1;
%     s=param(Nparam); % Swap rate
%     % Gaussian prior for the swap rate
%     p0(Nparam)=normpdf(s, 0.5, 1);
% end
% 
% % Construct joint distribution
% % Consider independent parameters here
% % We think it's generally acceptable for prior definition,
% % tho it's usually not the actual case
% p=1;
% for i=1:Nparam
%     p=p*p0(i);
% end
p=1;

end
