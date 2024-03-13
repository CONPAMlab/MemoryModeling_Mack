%% Group-Level Model Selection (SDs Model Recovery)
% Tianye Ma
% CONPAM lab (memory.ucr.edu)
% Department of Psychology, University of California, Riverside
%
%% Benchmark Models: Signal Detection

clc
close all
clear variables
%MA=MA2;
%----- Simulation -----%
disp("Now start simulation...")
Nsim=1;
Nsubj=100;
c=[-0.8 -0.4 0.4 0.8 1.2];
param_UVSD=[0.8, c, 1.1];
param_DPSD=[0.8, c, 1.1, 0.2, 0.05];
%param_DPSD=[MA.Param(1,:) 0.01];
% param_Mix=[1.5, c, 0.1, 0.8];
data_UVSD=cell(Nsim, length(Nsubj));
data_DPSD=cell(Nsim, length(Nsubj));
data_DPSD_perfect=cell(Nsim, length(Nsubj));
% data_Mix=cell(Nsim, length(Nsubj));
Ntrial=2000; % Ntrial per bin
% Ntarget=(length(c):-1:1)*Ntrial;
% Nlure=(length(c):-1:1)*Ntrial;
Ntarget = Ntrial/2;
Nlure = Ntrial/2;
for s=1:length(Nsubj)
    for i=1:Nsim
        i
        data0_all=cell(Nsubj(s),1);
        for subj=1:Nsubj(s)
            data0.Ntarget=Ntarget;
            data0.Nlure=Nlure;
            [data0.Nhit, data0.Nfa, r]=SimUVSD(param_UVSD, Ntarget, Nlure, 1, length(c));
            data0_all{subj}=data0;
        end
        data_UVSD{i,s}=data0_all;
    end
end
for s=1:length(Nsubj)
    for i=1:Nsim
        i
        data0_all=cell(Nsubj(s),1);
        for subj=1:Nsubj(s)
            data0.Ntarget=Ntarget;
            data0.Nlure=Nlure;
            [data0.Nhit, data0.Nfa, r]=SimDPSD(param_DPSD, Ntarget, Nlure, 1, length(c));
            data0_all{subj}=data0;
        end
        data_DPSD{i,s}=data0_all;
    end
end
% for s=1:length(Nsubj)
%     for i=1:Nsim
%         data0_all=cell(Nsubj(s),1);
%         for subj=1:Nsubj(s)
%             data0.Ntarget=Ntarget;
%             data0.Nlure=Nlure;
%             [data0.Nhit, data0.Nfa, r]=SimMixture(param_Mix, Ntarget, Nlure, 1, length(c));
%             data0_all{subj}=data0;
%         end
%         data_Mix{i,s}=data0_all;
%     end
% end
% % choose true model
% p_DPSD=0.7;
% data_all=cell(Nsim, length(Nsubj));
% for s=1:length(Nsubj)
%     for i=1:Nsim
%         N_DPSD=binornd(Nsubj(s), p_DPSD);  
%         
%     end
% end
disp("Done...")
disp("Now start building figures")
% UVSD
figure(1)
plot([0, 1], [0, 1], '--')
hold on
for i=1:Nsubj
    Nhit0=zeros(size(c));
    Nfa0=zeros(size(c));
    for j=1:Nsim
        Nhit0=Nhit0+data_UVSD{j}{i}.Nhit;
        Nfa0=Nfa0+data_UVSD{j}{i}.Nfa;
    end
    axis([0 1 0 1])
    plot(Nfa0/Nsim./Nlure, Nhit0./Ntarget/Nsim,'b*-')
end
% DPSD
figure(2)
plot([0, 1], [0, 1], '--')
hold on
for i=1:Nsubj
    Nhit0=zeros(size(c));
    Nfa0=zeros(size(c));
    for j=1:Nsim
        Nhit0=Nhit0+data_DPSD{j}{i}.Nhit;
        Nfa0=Nfa0+data_DPSD{j}{i}.Nfa;
    end
    axis([0 1 0 1])
    plot(Nfa0/Nsim./Nlure, Nhit0./Ntarget/Nsim,'b*-')
end
disp("Done...")

% % Create perfect datasets
% for s=1:length(Nsubj)
%     for i=1:Nsim
%         data0_all=cell(Nsubj(s),1);
%         for subj=1:Nsubj(s)
%             data0.Ntarget=Ntarget;
%             data0.Nlure=Nlure;
%             [data0.Nhit, data0.Nfa, r]=SimDPSD(param_DPSD, 1000*Ntarget, 1000*Nlure, 1, length(c));
%             data0.Nhit=round(data0.Nhit/1000);
%             data0.Nfa=round(data0.Nfa/1000);
%             data0_all{subj}=data0;
%         end
%         data_DPSD_perfect{i,s}=data0_all;
%     end
% end

%% Fit
% Names of the custom models should PERFECTLY match the function names
% FactorSpace={{'Unequal Variance'}, {'Unequal Variance', 'Oldness Threshold', 'Newness Threshold'}}; % UVSD & DPSD
FactorSpace={{'Unequal Variance','Recollection Threshold', 'Newness Threshold'}};
%FactorSpace={{'Unequal Variance'}};
Nmodel=length(FactorSpace);
MA_All=cell(1,Nmodel); % Pre-Allocation
for i=1:1
    fprintf('\n%s\n',FactorSpace{i}{1})
    Config_MA.Data=data_DPSD{1}; % e.g. E:/matlab/Data_BMW.mat
    % Specify model
    Config_MA.Model.Model='SignalDetection_ROC';
    % Model Variants
    Config_MA.Model.Variants=FactorSpace{i}; % loop factors
    % Specity mode
    Config_MA.FitOptions.Method='MLE';
    
%     % The initial values and the constraints of the custom models need to be manually assigned
%     % UVSD
%     c_start=[-0.3 -0.2 0.1 0.2 0.3];
%     Config_MA.Constraints.start=[1 c_start 1]; % Initial guess
%     Config_MA.Constraints.lb=[-0.5, -3*ones(1,length(c_start)), 0.0001];
%     Config_MA.Constraints.ub=[10, 3*ones(1,length(c_start)), 10];
%     
    % The initial values and the constraints of the custom models need to be manually assigned
    % DPSD
    c_start=[-0.3 -0.2 0.1 0.2 0.3];
    Config_MA.Constraints.start=[1 c_start 1 0.2 0.1]; % Initial guess
    Config_MA.Constraints.lb=[-0.5, -3*ones(1,length(c_start)), 0.0001, 0.0001 0.0001];
    Config_MA.Constraints.ub=[10, 3*ones(1,length(c_start)),10 1 1];
    
    Config_MA.FitOptions.Algorithm='fmincon: sqp'; % Optimization algorithm
    % Configuration
    MA=Configuration_BMW(Config_MA);
    % Estimation
    MA=ModelFit_BMW(MA);
    MA_All{i}=MA;
end

%% Check Goodness-of-Fit
figure(3)
for data_id=1:2
hold on
plot([0, 1], [0, 1], '--')

% Nhit1=Config_MA.Data{data_id}.Nhit;
% Nfa1=Config_MA.Data{data_id}.Nfa;
% axis([0 1 0 1])
% plot(Nfa1/Ntrial, Nhit1/Ntrial,'r*-')

% recovery
Nsim=100;
Nsubj=10;
param_check=MA.Param(data_id,:);
data_check=cell(Nsim, length(Nsubj));
Ntarget=Ntrial*ones(length(c),1);
Nlure=Ntrial*ones(length(c),1);

for s=1:length(Nsubj)
    for i=1:Nsim
        data0_all=cell(Nsubj(s),1);
        for subj=1:Nsubj(s)
            data0.Ntarget=Ntarget;
            data0.Nlure=Nlure;
            [data0.Nhit, data0.Nfa, r]=SimDPSD(param_check, Ntarget, Nlure, 1, length(c)); % change the function name
            data0_all{subj}=data0;
        end
        data_check{i,s}=data0_all;
    end
end

pfa_all=zeros(Nsubj,length(c));
phit_all=zeros(Nsubj,length(c));
for i=1:Nsubj
    Nhit0=zeros(size(c));
    Nfa0=zeros(size(c));
    for j=1:Nsim
        Nhit0=Nhit0+data_check{j}{i}.Nhit;
        Nfa0=Nfa0+data_check{j}{i}.Nfa;
    end
    axis([0 1 0 1])
    pfa_all(i,:)=Nfa0/Nsim/Ntrial;
    phit_all(i,:)=Nhit0/Ntrial/Nsim;
end
plot(mean(pfa_all), mean(phit_all),'b*-')

end


