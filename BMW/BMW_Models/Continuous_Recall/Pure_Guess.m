%% Pure Guess (van den Berg et al., 2012)
%
% Define the Slots-plus-Averaging model
% ------------ 
% Output=Slots_plus_Averaging(param, Data, Input)
%
% ## Theory ##
% This model assumed that there are a limited number of resouce slots and
% the fidelity of each representation depends on the number of slots being
% allocated to the corresponding item.
% The response probability density distribution is the weighted addition
% between random guess (uniform) & memory response (Von Mises)
%
% ## Input ##
% check the manual for details (BMW('manual'))
%
% - param
% K, kappa_1, kappa_r, (bias, biasF, precF, s)
%
% - Data
% Data.error (response-sample), Data.SS (set size)
% (Data.sample_range, Data.sample, Data.error_nt, Data.sample_nt)
%
% - Input
% Input.Variant
%   options of model variants
%       Input.Variants.Bias, 0/1 to decide whether consider representational
%       shift/response bias
%       Input.Variants.Swap, 0/1 to decide whether use swap variants
%       Input.Variants.BiasF, 0/1 to decide whether consider a cosine-shaped
%       fluctuation of representational shift/response bias
% Input.Output
%   string, choose output mode
%       'Prior', only output prior density
%       'LLH', output log likelihood
%       'LP', output log posterior density
% Input.PDF
%   0/1, output pdf or not. default as 0
%   Valid only when Input.Output=='LLH'
%
% ## Output ##
% Output is conditional to Input.Output & Input.PDF
%
% ## Reference ##
% - Zhang, W., & Luck, S. J. (2008). "Discrete fixed-resolution representations in visual working memory".
% Nature, 453(7192), 233.
% - van den Berg, R., Shin, H., Chou, W. C., George, R., & Ma, W. J. (2012).
% "Variability in encoding precision accounts for visual short-term memory limitations".
% Proceedings of the National Academy of Sciences, 109(22), 8780-8785.
% - Bays, P. M., Catalao, R. F., & Husain, M. (2009). "The precision of visual working memory
% is set by allocation of a shared resource." Journal of Vision, 9(10), 7-7.
% - Pratte, M. S., Park, Y. E., Rademaker, R. L., & Tong, F. (2017). "Accounting for stimulus-specific variation
% in precision reveals a discrete capacity limit in visual working memory."
% Journal of Experimental Psychology: Human Perception and Performance, 43(1), 6.
%
% ------------
% Programmed by Ma, Tianye
% Under the guidance of Dr. Ku, Yixuan
% Memory, Attention & Cognition (MAC) Lab,
% East China Normal University
% 9/26/2019
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function Output=Pure_Guess(param, Data, Input)

    % Specify parameters
    Nparam=0;
    if ~isfield(Input,'Variants') % No Variants
        Input.Variants={};
    end
%     if any(strcmp(Input.Variants,'ResponseNoise'))
%         Nparam=Nparam+1;
%         kappa_r=param(Nparam); % Response precision
%     end
    
    % Configuration
    samples=Data.sample;
    responses=Data.response;
    error_range=Data.error_range;
    if length(error_range)==2
        continuous=1;
        period=error_range(2)-error_range(1);
    else
        continuous=0;
        period=max(error_range)-min(error_range)+(error_range(2)-error_range(1));
    end
    errors=CircDist_BMW('Diff',responses,samples,period);
    SS=Data.SS;
    SS_range=unique(SS);
    
    if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
        Prior=prior(param, Input); % get prior
    elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
        Prior=1; % uniform prior
    end
    
    if ~strcmp(Input.Output,'Prior')
        % LH function
        if continuous==1
            p_error=zeros(1,length(errors));
            p_error_NT=zeros(1,length(errors));
            for i_error=1:length(errors)
                        p_error(i_error)=1/(2*period/360)/pi;    
            end
            p_T=(1-s)*p_error;
            p_NT=s*p_error_NT;
            p_LH=p_T+p_NT;
        else

            p_error=zeros(length(SS_range), length(error_range),length(sample_range));
            for i_N=1:length(SS_range)
                for i_error=1:length(error_range)
                    p_error(i_N,i_error,:)=1/(2*period/360)/pi;
                end
                for i_sample=1:length(sample_range)
                    p_error(i_N,:,i_sample)=p_error(i_N,:,i_sample)/sum(p_error(i_N,:,i_sample));
                end
            end
            
            % Calculate LH
            p_T=zeros(1,length(errors));
            p_NT=zeros(1,length(errors));
            for i=1:length(errors)
                    p_T(i)=(1-s)*p_error(SS_range==SS(i),error_range==errors(i), 1);
            end
            p_LH=p_T+p_NT; % Target + non-target LH
        end
        
        % LLH
        if isfield(Input,'PDF') && Input.PDF==1
            LLH.error=p_error; % PDF
        else
            LLH=-sum(log(p_LH)); % Negative LLH
        end
        
        % Posterior
        LP=-log(Prior)+LLH; % likelihood*prior
        
        if LP==Inf || isnan(LP)
            LP=realmax('double'); % Output should be a real value
        end
        
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

end

% Define prior
function p=prior(param, Input)

% Specify parameters
K=param(1); % Capacity
% weibull prior for capacity
p0(1)=wblpdf(K,3.5,3); % given that K is ofter 3~4
kappa_1=param(2); % Unit resource
% Gamma prior for unit resource
p0(2)=gampdf(kappa_1,3,5);
Nparam=2;

if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
if any(strcmp(Input.Variants,'ResponseNoise'))
    kappa_r=param(3); % Response variability
    % Gamma prior for response noise
    p0(3)=gampdf(kappa_r,3,5);
end
if any(strcmp(Input.Variants,'Bias'))
    Nparam=Nparam+1;
    bias=param(Nparam); % Mean bias
    % Gaussian prior for bias
    p0(Nparam)=normpdf(bias, 0, 1);
end
if any(strcmp(Input.Variants,'BiasF'))
    Nparam=Nparam+1;
    biasF=param(Nparam); % Fluctuation of bias
    % Gaussian prior for the fluctuation of bias
    p0(Nparam)=normpdf(biasF, 0, 5);
end
if any(strcmp(Input.Variants,'PrecF'))
    Nparam=Nparam+1;
    precF=param(Nparam); % Fluctuation of precision
    % Gaussian prior for the fluctuation of precision
    p0(Nparam)=normpdf(precF, 0, 1);
end
if any(strcmp(Input.Variants,'Swap'))
    Nparam=Nparam+1;
    s=param(Nparam); % Swap rate
    % Gaussian prior for the swap rate
    p0(Nparam)=normpdf(s, 0.5, 1);
end

% Construct joint distribution
% Consider independent parameters here
% We think it's generally acceptable for prior definition,
% tho it's usually not the actual case
p=1;
for i=1:Nparam
    p=p*p0(i);
end

end
