% SimMixture(param, Ntarget, Nlure, Nsubj, Nc)
% SimDPSD(param, Ntarget, Nlure, Nsubj, Nc)

dp_range=[1 3 5];
c=[-1.0 -0.4 0.4 1.0 1.8];
Ro_range=[0.01 0.1:0.:0.9];
sigma=1;
Rn = 0;
Nsubj = 20;
Ntarget = 2000;
Nlure = 2000;

Data_DPSD = cell(length(dp_range), length(Ro_range));
for dp_id = 1:length(dp_range)
    for Ro_id = 1:length(Ro_range)
        param_DPSD = [dp_range(dp_id), c, sigma, Ro_range(Ro_id), Rn];
        [Data_DPSD{dp_id, Ro_id}] = SimDPSD(param_DPSD, Ntarget, Nlure, Nsubj, length(c));
    end
end

% Create a list of model factors
% allFactor = {'Unequal Variance', 'Recollection Threshold', 'Newness Threshold', 'Recognition Rate'}; 
allFactor = {'Recognition Rate'}; 
% allFactor = {'Recollection Threshold', 'Unequal Variance'}; 
% allFactor = {'Recollection Threshold'}; 
Nfactor = length(allFactor);
all_fc = de2bi(1:(2^Nfactor));
all_fc = all_fc(:, 1:Nfactor);
all_fc = all_fc(1:2,:);
Factor_List= cell(size(all_fc, 1), 1);

for f = 1:size(all_fc, 1)
    Factor_List{f} = allFactor(logical(all_fc(f,:)));
end
% The same model: SDT
Model_List = repmat({'Signal Detection (ROC)'}, size(all_fc, 1), 1);
% Fit
FitResults = cell(length(dp_range), length(Ro_range));
for dp_id = 1:length(dp_range)
    for Ro_id = 1:length(Ro_range)
        Config_MA.Data = Data_DPSD{dp_id, Ro_id};
        FitResults_c = cell(size(all_fc, 1), 1);
        for f_id = 1:size(all_fc, 1)
            
            if isempty(cell2mat(Factor_List{f_id}))
                Variants_Display='None';
            else
                Variants_Now=Factor_List{f_id};
                Variants_Display=cell(1,length(Variants_Now));
                for j=1:length(Variants_Now)
                    Variants_Display{j}=[Variants_Now{j},' '];
                end
                Variants_Display=cell2mat(Variants_Display);
            end
            
            fprintf('\ndp: %d, Ro: %d, Factor: %s\n', dp_id, Ro_id, Variants_Display)
            Config_MA.Model.Model = Model_List{f_id};
            Config_MA.Model.Variants = Factor_List{f_id};
            
            Config_MA.FitOptions.Method='MLE';
            Config_MA.FitOptions.Algorithm='fmincon: sqp'; % Optimization algorithm
            Config_MA.Criteria={'BIC','AIC','AICc','LLH'};
            Config_MA.Model.Ntrial = Ntarget;
%             % The initial values and the constraints of the custom models need to be manually assigned
%             c_start=[-0.3 -0.2 0.1 0.2 0.3];
%             Config_MA.Constraints.start=[1 c_start 1]; % Initial guess
%             Config_MA.Constraints.lb=[-0.5, -3*ones(1,length(c_start)), 0.0001];
%             Config_MA.Constraints.ub=[10, 3*ones(1,length(c_start)), 10];
%             
            MA=Configuration_BMW(Config_MA);
            FitResults_c{f_id}=ModelFit_BMW(MA); % Run!

        end
            FitResults{dp_id, Ro_id} = FitResults_c;
    end
end

Pm_now = zeros(length(dp_range), length(Ro_range));
Ro_now = zeros(length(dp_range), length(Ro_range));
dpSlot_now = zeros(length(dp_range), length(Ro_range));
for i = 1:length(dp_range)
    for j = 1:length(Ro_range)
        Pm_now(i,j) = median(FitResults{i,j}{1}.Param(:,end));
%         Ro_now(i,j) = median(FitResults{i,j}{1}.Param(:,end));
        dpSlot_now(i,j) = median(FitResults{i,j}{1}.Param(:,1));
    end
end