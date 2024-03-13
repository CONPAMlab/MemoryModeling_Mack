%% Target Confusability Competition
%
%
% Define the likelihood function for the TCC model
% ------------
% Output=TCC(param, Data, Input)
%
% ## Theory ##
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
%
% ------------
% Programmed by Ma, Tianye
% 10/23/2023
%
% Bug reports or any other feedbacks please contact M.T. (mack_ma2018@outlook.com)
% BMW toolbox: https://github.com/Mack-Ma/Bayesian_Modeling_of_Working_Memory
%

function Output=TCC(param, Data, Input)
% Specify parameters
SS=Data.SS;
SS_range=unique(SS);
gain = param(1); % memory strength
dprime = gain./SS_range;
Nparam = 1;
if ~isfield(Input,'Variants') % No Variants
    Input.Variants={};
end
if any(strcmp(Input.Variants,'Precision'))
    Nparam=Nparam+1;
    kappa = param(Nparam);
    confusionVector = 0;
    FreeCV = 1;
else
    confusionVector = Data.ConfusionVector;
    FreeCV = 0;
    kappa = 1;
end
if any(strcmp(Input.Variants,'Capacity'))
    Nparam=Nparam+1;
    K=param(Nparam); 
    g = 1-K./SS_range;
    g(g<0) = 0;
    g(g>1) = 1;
else
    g=zeros(length(SS_range),1);
end

% Configuration
samples=Data.sample;
responses=Data.response;
error_range=Data.error_range;
period=max(error_range)-min(error_range)+(error_range(2)-error_range(1));
errors=CircDist_BMW('Diff',responses,samples,period);

if strcmp(Input.Output,'LP') || strcmp(Input.Output,'Prior') || strcmp(Input.Output,'All')
    Prior=prior(param, Input); % get prior
elseif strcmp(Input.Output,'LLH') || strcmp(Input.Output,'LPPD')
    Prior=1; % uniform prior
end


if ~strcmp(Input.Output,'Prior')
        
    s_range = -10:0.05:10;
    confusionVector = confusionVector + FreeCV*exp(-kappa.*(360/period)*abs(error_range));
% Using the psychological distance data, make a CDF for each possible
    % error amount (-180 to 180) in memory strength space 
    maxDist = ones(size(s_range));
    for i=1:length(error_range)
      % How likely is an item of every possible psych. distance to generate
      % every possible memory strength signal:
      normals(i,:) = normcdf(s_range,...
        confusionVector(i).*dprime,1);
      
      % How likely is each psych. distance to generate the maximum memory
      % match signal? Distribution of the maximums is just prod(CDFs):
      maxDist = maxDist .* normals(i,:);
    end

    pdf = zeros(length(error_range),1);
    for j=1:period
      % Remove the current signal from the max distribution(we want to know
      % the chance that ANOTHER color was larger than this one) and turn
      % this distribution into a PDF instead of a cdf:
      maxDistCur = maxDist ./ normals(j,:);
      maxPdf = [0 diff(maxDistCur)];
      maxPdf = maxPdf ./ nansum(maxPdf);
      
      % Now ask the chance this value was greater than all of the other
      % ones:
      chanceOverThisVal = (1-normals(j,:));
      pdf(j) = nansum(maxPdf.*chanceOverThisVal);
    end   
    
    motorError = 2;
    widthToConsider = motorError * 3;
    gaussFilter = normpdf(-widthToConsider:widthToConsider, 0, 2);
    gaussFilter = gaussFilter / sum(gaussFilter);
    
    % Take this pdf and put it three times, then convolve it with motor
    % noise, than pull out the center part (so as to avoid trouble with
    % wrapping around the circle in the motor noise part):
    c = [pdf' pdf' pdf'];
    c = conv(c, gaussFilter, 'same');
    pdf = c((length(error_range)+1):(2*length(error_range)))';
   
    % Now normalize resulting distribution:
    pdf = pdf./sum(pdf);
    pdf = pdf .* (1-g) + 1/360*g;
    pdf = pdf ./ sum(pdf);
    
    % Now fill in the probability for each error in the data by replacing
    % it with the relevant amount from the PDF:
    c = floor(errors+length(error_range)/2);
    c(c==0) = period;
    p = nan(size(errors));
    p(~isnan(c)) = pdf(c(~isnan(c)));
    p_LH = p;
    
    % LLH
    if isfield(Input,'PDF') && Input.PDF==1
        LLH.error=p_error; % PDF
    else
        LLH=-sum(log(p_LH)); % Negative LLH
    end

    if LLH==Inf || isnan(LLH)
        LLH=realmax('double'); % Output should be a real value
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


function Prior = prior(param, Input)
    Prior = 1;
end
% function p=strength_pdf(s, mu1, mu2)
% p=normpdf(s, mu1, 1); % pick one strength
% for i=1:length(mu2)
%     p_max=normcdf(s, mu2(i), 1); % probability of being larger than the other options
%     p=p.*p_max;
% end
% end