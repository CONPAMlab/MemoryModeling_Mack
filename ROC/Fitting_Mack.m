%% Fitting SDT models
% TM@UCR
% Dec 2023

%% Load and merge data
clc
clear
close all

addpath('/home/tianye/Matlab_WS/Ensemble_Inik')

SDT_dlist = {'DataAll_ROC_Pk_Exp3', 'DataAll_ROC_New'};
% SDT_dlist = {'DataAll_ROC_New'};
Ndata = length(SDT_dlist);
SDT_AllC = cell(Ndata, 1);
cond_list = {[1 2], [1 2]}; % Specify task condition
Nrating = 6;

for i = 1:Ndata
    dname = SDT_dlist{i};
    cond0 = cond_list{i};
    load(dname)
    Nsubj = length(DataAll);
    datac = cell(length(cond0), 1);
    for c = 1:length(cond0)
        data0 = cell(Nsubj, 1);
        for subj = 1:Nsubj
            cond_now = cond0(c);
            for r = 1:Nrating
                data_haha.Nhit(r) = sum(DataAll{subj}.IsSame(DataAll{subj}.Rating>=r & DataAll{subj}.Task == cond_now)==0);
                data_haha.Nfa(r) = sum(DataAll{subj}.IsSame(DataAll{subj}.Rating>=r & DataAll{subj}.Task == cond_now)==1);
            end
            data0{subj}=data_haha;
        end
        datac{c} = data0;
    end
    SDT_AllC{i} = datac;
end

%% Fit the models
% Create a list of model factors
% allFactor = {'Unequal Variance', 'Recollection Threshold', 'Newness Threshold', 'Recognition Rate'}; 
% allFactor = {'Unequal Variance', 'Recollection Threshold', 'Newness Threshold'}; 
% allFactor = {'Recollection Threshold', 'Unequal Variance'}; 
allFactor = {'Recollection Threshold'}; 
Nfactor = length(allFactor);
all_fc = de2bi(1:(2^Nfactor));
all_fc = all_fc(:, 1:Nfactor);
Factor_List= cell(size(all_fc, 1), 1);
for f = 1:size(all_fc, 1)
    Factor_List{f} = allFactor(logical(all_fc(f,:)));
end
% The same model: SDT
Model_List = repmat({'Signal Detection (ROC)'}, size(all_fc, 1), 1);

% Fit
FitResults = cell(Ndata, 1);
for d_id = 1:Ndata
    cond0 = cond_list{d_id};
    FitResults_d = cell(length(cond0), 1);
    for c_id = 1:length(cond0)
        Config_MA.Data = SDT_AllC{d_id}{c_id};
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
            
            fprintf('\nData: %d, Condition: %d, Factor: %s\n', d_id, c_id, Variants_Display)
            Config_MA.Model.Model = Model_List{f_id};
            Config_MA.Model.Variants = Factor_List{f_id};
            
            Config_MA.FitOptions.Method='MLE';
            Config_MA.FitOptions.Algorithm='fmincon: sqp'; % Optimization algorithm
            Config_MA.Criteria={'BIC','AIC','AICc','LLH'};
            Config_MA.Model.Ntrial = 150;
%             % The initial values and the constraints of the custom models need to be manually assigned
%             c_start=[-0.3 -0.2 0.1 0.2 0.3];
%             Config_MA.Constraints.start=[1 c_start 1]; % Initial guess
%             Config_MA.Constraints.lb=[-0.5, -3*ones(1,length(c_start)), 0.0001];
%             Config_MA.Constraints.ub=[10, 3*ones(1,length(c_start)), 10];
%             
            MA=Configuration_BMW(Config_MA);
            FitResults_c{f_id}=ModelFit_BMW(MA); % Run!

        end
        FitResults_d{c_id} = FitResults_c;
    end
    FitResults{d_id} = FitResults_d;
end

%% Model Comparison
% Group-level model comparison
MC = cell(Ndata, 1);
for d_id = 1:Ndata
    cond0 = cond_list{d_id};
    MC_d = cell(length(cond0), 1);
    for c_id = 1:length(cond0)
        MC_d{c_id} =  ModelComparison_BMW(FitResults{d_id}{c_id}');
    end
    MC{d_id} = MC_d;
end
% Family-based model comparison
MC_fam = cell(Ndata, 1);
contrast = all_fc(:,1) + 1; %
for d_id = 1:Ndata
    cond0 = cond_list{d_id};
    MC_d = cell(length(cond0), 1);
    for c_id = 1:length(cond0)
        MC_d{c_id} =  ModelComparison_BMW(FitResults{d_id}{c_id}', contrast);
    end
    MC_fam{d_id} = MC_d;
end

%% Visualization
% figure(1)
% for d_id = 1:Ndata
%     subplot(1,2,1)
% end


%%
rmpath('/home/tianye/Matlab_WS/Ensemble_Inik')

