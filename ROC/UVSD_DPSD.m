%% Test hahaha
% Simply UVSD vs. DPSD (without Rn)

%% Prologue

close all
clear
clc

%----- Simulation -----%
disp("Now start simulation...")
% Nsim=100;
Nsubj_range=[25];
c=[-0.8 -0.4 0.4 0.8 1.2];
Ro_range=[0.01 0.3 0.7];
dp_range=[0.001 1.5 3.5];
sigma=1.1;
Ntrial=200;
Ntarget=(length(c):-1:1)*Ntrial;
Nlure=(length(c):-1:1)*Ntrial;

BIC_1=cell(length(Nsubj_range), length(Ro_range));
BIC_2=cell(length(Nsubj_range), length(Ro_range));
data1=cell(length(Nsubj_range), length(Ro_range));
data2=cell(length(Nsubj_range), length(Ro_range));
r_BMS=zeros(length(Nsubj_range), length(Ro_range));
r_PMP=zeros(length(Nsubj_range), length(Ro_range));
r_per=zeros(length(Nsubj_range), length(Ro_range));
scaling=1000;
for s=1:length(dp_range)
    for ro=1:length(Ro_range)
        
        % simulation
        p_1=0.6; % p_UVSD
        Nsubj=Nsubj_range(1);
        % simulate UVSD
        data1_all=cell(ceil(Nsubj*p_1),1);
        for subj=1:ceil(Nsubj*p_1)
            data0.Ntarget=Ntarget;
            data0.Nlure=Nlure;
            [data0.Nhit, data0.Nfa, r]=SimUVSD([dp_range(s), c, sigma], Ntarget*scaling, Nlure*scaling, 1, length(c));
            data0.Nhit=round(data0.Nhit/scaling);
            data0.Nfa=round(data0.Nfa/scaling);
            data1_all{subj}=data0;
        end
        data1{s,ro}=data1_all;
        % simulate DPSD
        data2_all=cell(Nsubj-ceil(Nsubj*p_1),1);
        for subj=1:Nsubj-ceil(Nsubj*p_1)
            data0.Ntarget=Ntarget;
            data0.Nlure=Nlure;
            [data0.Nhit, data0.Nfa, r]=SimDPSD([dp_range(s), c, sigma, Ro_range(ro), 0], Ntarget*scaling, Nlure*scaling, 1, length(c));
            data0.Nhit=round(data0.Nhit/scaling);
            data0.Nfa=round(data0.Nfa/scaling);
            data2_all{subj}=data0;
        end
        data2{s,ro}=data2_all;
        % merge
        data_all=[data1_all; data2_all];
        
        % fit
        Factor_UVSD={'Unequal Variance'};
        Factor_DPSD={'Recollection Threshold'};
        Config_MA.Model.Model='SignalDetection_ROC'; % Specify model
        Config_MA.FitOptions.Method='MLE';
        Config_MA.FitOptions.Algorithm='fmincon: sqp'; % Optimization algorithm
        Config_MA.Criteria={'BIC','AIC','AICc','LLH'};
        Config_MA.Ntrial=Ntrial;
        % fit UVSD
        Config_MA.Data=data_all; % e.g. E:/matlab/Data_BMW.mat
        Config_MA.Model.Variants=Factor_UVSD; % factors
        % The initial values and the constraints of the custom models need to be manually assigned
        c_start=[-0.3 -0.2 0.1 0.2 0.3];
        Config_MA.Constraints.start=[1 c_start 1]; % Initial guess
        Config_MA.Constraints.lb=[-0.5, -3*ones(1,length(c_start)), 0.0001];
        Config_MA.Constraints.ub=[10, 3*ones(1,length(c_start)), 10];
        MA1_0=Configuration_BMW(Config_MA);
        MA1=ModelFit_BMW(MA1_0); % run!
        % fit DPSD
        Config_MA.Data=data_all; % e.g. E:/matlab/Data_BMW.mat
        Config_MA.Model.Variants=Factor_DPSD; % factors
        % The initial values and the constraints of the custom models need to be manually assigned
        c_start=[-0.3 -0.2 0.1 0.2 0.3];
        Config_MA.Constraints.start=[1 c_start 0.2]; % Initial guess
        Config_MA.Constraints.lb=[-0.5, -3*ones(1,length(c_start)), 0.0001];
        Config_MA.Constraints.ub=[10, 3*ones(1,length(c_start)) 1];
        MA2_0=Configuration_BMW(Config_MA);
        % Estimation
        MA2=ModelFit_BMW(MA2_0);
        % get BICs
        BIC_1{s,ro}=MA1.AIC;
        BIC_2{s,ro}=MA2.AIC;
        
        % Group-level comparison
        BIC_all=[BIC_1{s,ro},BIC_2{s,ro}]';
        BIC_all=BIC_all-repmat(min(BIC_all,[],1),[2,1]);
        % RFX-BMS
        Opt_BMC.Start=1e-6*ones(2,1); % Flat prior
        Opt_BMC.MaxIter=1e6;
        Opt_BMC.Stop=1e-6;
        Opt_BMC.Verbosity=0;
        Opt_BMC.Rec=0;
        BMS_all=BMW_BMS(-BIC_all, Opt_BMC);
        r_BMS(s,ro)=BMS_all.r(1);
        % mean posterior model probability p(model | data)
        pmp_1=exp(-.5*BIC_all(1,:))./sum(exp(-.5*BIC_all));
        pmp_2=exp(-.5*BIC_all(2,:))./sum(exp(-.5*BIC_all));
        r_PMP(s,ro)=mean(pmp_1);
        % percentage of being the best model
        r_per(s,ro)=mean(BIC_1{s,ro}<=BIC_2{s,ro});
    end
end

figure(1)
subplot(3,1,1)
plot(r_BMS'-p_1, 'b.-')
hold on
plot([1, length(dp_range)],[0 0])
subplot(3,1,2)
plot(r_PMP'-p_1, 'b.-')
hold on
plot([1, length(dp_range)],[0 0])
subplot(3,1,3)
plot(r_per'-p_1, 'b.-')
hold on
plot([1, length(dp_range)],[0 0])

figure(2)
for ro=1:length(Ro_range)
    for s=1:length(Nsubj_range)
        LW=s*0.5;
        color_line=[0.3*ro, 0, 0];
        plot(BIC_1{s,ro}-BIC_2{s,ro},'Color',color_line,'LineWidth',LW)
        hold on
    end
end
plot([0, max(Nsubj_range)],[0,0],'--','Color',[0 0 0])
