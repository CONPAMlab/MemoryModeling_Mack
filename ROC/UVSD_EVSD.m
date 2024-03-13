%% Test hahaha
% Simply UVSD vs. EVSD

%% Prologue

close all
clear
clc

%----- Simulation -----%
disp("Now start simulation...")
% Nsim=100;
Ntrial_range=100:100:800;
c=[-0.9 -0.4 0 0.4 0.9];
Vo_range=0.5:0.5:4;
dp=1.5;
p_range=[0.7 0.3];
% sigma=1.1;
Nsubj=50;
Ntarget=(length(c):-1:1)'*Ntrial_range/length(c);DataHAHA = cell(length(subj_id));

Nlure=(length(c):-1:1)'*Ntrial_range/length(c);

BIC_1=cell(length(Ntrial_range), length(Vo_range));
BIC_2=cell(length(Ntrial_range), length(Vo_range));
data1=cell(length(Ntrial_range), length(Vo_range));
data2=cell(length(Ntrial_range), length(Vo_range));
r_BMS=zeros(length(Ntrial_range),length(Vo_range),length(p_range));
r_ME=zeros(length(Ntrial_range),length(Vo_range),length(p_range));
r_ME2=zeros(length(Ntrial_range),length(Vo_range),length(p_range));
r_PMP=zeros(length(Ntrial_range),length(Vo_range),length(p_range));
r_PMP2=zeros(length(Ntrial_range),length(Vo_range),length(p_range));
r_per=zeros(length(Ntrial_range),length(Vo_range),length(p_range));
scaling=1;
for haha=1:length(p_range)
    for s=1:length(Ntrial_range)
        for ro=1:length(Vo_range)
            
            % simulation
            p_1=p_range(haha); % p_UVSD
            % simulate UVSD
            data1_all=cell(ceil(Nsubj*p_1),1);
            for subj=1:ceil(Nsubj*p_1)
                data0.Ntarget=Ntarget(:,s);
                data0.Nlure=Nlure(:,s);
                [data0.Nhit, data0.Nfa, r]=SimUVSD([dp, c, Vo_range(ro)], Ntarget(:,s)*scaling, Nlure(:,s)*scaling, 1, length(c));
                data0.Nhit=round(data0.Nhit/scaling);
                data0.Nfa=round(data0.Nfa/scaling);
                data1_all{subj}=data0;
            end
            data1{s,ro}=data1_all;
            % simulate EVSD
            data2_all=cell(Nsubj-ceil(Nsubj*p_1),1);
            for subj=1:Nsubj-ceil(Nsubj*p_1)
                data0.Ntarget=Ntarget(:,s);
                data0.Nlure=Nlure(:,s);
                [data0.Nhit, data0.Nfa, r]=SimUVSD([dp, c, 1], Ntarget(:,s)*scaling, Nlure(:,s)*scaling, 1, length(c));
                data0.Nhit=round(data0.Nhit/scaling);
                data0.Nfa=round(data0.Nfa/scaling);
                data2_all{subj}=data0;
            end
            data2{s,ro}=data2_all;
            % merge
            data_all=[data1_all; data2_all];
            
            % fit
            Factor_UVSD={'Unequal Variance'};
            Factor_EVSD={};
            Config_MA.Model.Model='SignalDetection_ROC'; % Specify model
            Config_MA.FitOptions.Method='MLE';
            Config_MA.FitOptions.Algorithm='fmincon: sqp'; % Optimization algorithm
            Config_MA.Criteria={'BIC','AIC','AICc','LLH'};
            Config_MA.Ntrial=Ntrial_range(s);
            % fit UVSD
            Config_MA.Data=data_all; % e.g. E:/matlab/Data_BMW.mat
            Config_MA.Model.Variants=Factor_UVSD; % factors
            % The initial values and the constraints of the custom models need to be manually assigned
            c_start=[-0.4 -0.2 0 0.2 0.4];
            Config_MA.Constraints.start=[1 c_start 1]; % Initial guess
            Config_MA.Constraints.lb=[-0.5, -3*ones(1,length(c_start)), 0.0001];
            Config_MA.Constraints.ub=[10, 3*ones(1,length(c_start)), 10];
            MA1_0=Configuration_BMW(Config_MA);
            MA1=ModelFit_BMW(MA1_0); % run!
            % fit EVSD
            Config_MA.Data=data_all; % e.g. E:/matlab/Data_BMW.mat
            Config_MA.Model.Variants=Factor_EVSD; % factors
            % The initial values and the constraints of the custom models need to be manually assigned
            c_start=[-0.4 -0.2 0 0.2 0.4];
            Config_MA.Constraints.start=[1 c_start]; % Initial guess
            Config_MA.Constraints.lb=[-0.5, -3*ones(1,length(c_start))];
            Config_MA.Constraints.ub=[10, 3*ones(1,length(c_start))];
            MA2_0=Configuration_BMW(Config_MA);
            % Estimation
            MA2=ModelFit_BMW(MA2_0);
            % get BICs
            BIC_1{s,ro}=MA1.BIC;
            BIC_2{s,ro}=MA2.BIC;
            
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
            r_BMS(s,ro,haha)=BMS_all.r(1);
            % mean/median model evidence
            BIC_mean=mean(BIC_all,2);
            BIC_mean=BIC_mean-repmat(min(BIC_mean,[],1),[2,1]);
            BIC_mean2=log(mean(exp(BIC_all),2));
            r_ME(s,ro,haha)=exp(-.5*BIC_mean(1,:))./sum(exp(-.5*BIC_mean));
            r_ME2(s,ro,haha)=exp(-.5*BIC_mean2(1,:))./sum(exp(-.5*BIC_mean2));
            % mean/median posterior model probability, p(model | data)
            pmp_1=exp(-.5*BIC_all(1,:))./sum(exp(-.5*BIC_all));
            pmp_2=exp(-.5*BIC_all(2,:))./sum(exp(-.5*BIC_all));
            r_PMP(s,ro,haha)=mean(pmp_1);
            r_PMP2(s,ro,haha)=median(pmp_1);
            % percentage of being the best model
            r_per(s,ro,haha)=mean(BIC_1{s,ro}<=BIC_2{s,ro});
        end
    end
end

% Visualization
% colorM = [0.01*ones(length(Ntrial_range),1),(0.1:(0.9/length(Ntrial_range)):(1-0.9/length(Ntrial_range)))',  0.01*ones(length(Ntrial_range),1)];
colorM0=colormap('autumn'); % Get color matrix
color_id=floor(linspace(1,256,length(Ntrial_range)));
colorM=colorM0(color_id,:);
hi=2;
fighaha=figure(1);
set(fighaha,'Position',[71,73,1500,756])
for hei=1:length(Ntrial_range)
    %---%
    subplot(3,2,6)
    title('RFX-BMS')
    if all(abs(r_BMS(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 length(Vo_range)+1])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticks('')
    plot(r_BMS(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(Vo_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,2)
    title('Mean PMP')
    if all(abs(r_PMP(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 length(Vo_range)+1])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticks('')
    plot(r_PMP(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(Vo_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,4)
    title('Median PMP')
    if all(abs(r_PMP2(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 length(Vo_range)+1])
    xticks('')
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticklabels({num2str(Vo_range(1)), num2str(Vo_range(round(length(Vo_range)/2))), num2str(Vo_range(end))})
    plot(r_PMP2(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(Vo_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,5)
    title('Best Model Percentage')
    xlim([0 length(Vo_range)+1])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xlabel('Vo')
    xticks([1, length(Vo_range)])
    xticklabels({'0.5','4.0'})
    if all(abs(r_per(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    plot(r_per(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(Vo_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,1)
    title('Mean Log Model Evidence')
    ylabel('Error of Recovery')
    if all(abs(r_ME(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 length(Vo_range)+1])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticks('')
    plot(r_ME(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(Vo_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,3)
    title('Mean Model Evidence')
    if all(abs(r_ME2(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 length(Vo_range)+1])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticks('')
    plot(r_ME2(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(Vo_range)],[0 0],'Color',[0,0,0])
end
cb=colorbar('Location','layout',...
    'Position',[0.9247,0.1111,0.0153,0.8108],...
    'Xtick',[0, 0.5, 1],...
    'XtickLabel',{'10', '100', '200'});
cb.Label.String='Sample Size';
cb.Label.FontSize=11;

%
hi=1
err{5}=r_per(:,:,hi)'-p_range(hi);
err{6}=r_BMS(:,:,hi)'-p_range(hi);
err{3}=r_PMP(:,:,hi)'-p_range(hi);
err{4}=r_PMP2(:,:,hi)'-p_range(hi);
err{1}=r_ME(:,:,hi)'-p_range(hi);
err{2}=r_ME2(:,:,hi)'-p_range(hi);
errM1=zeros(6,6);
for haha=1:6
    for hihi=1:6
        errM1(haha,hihi)=sum(sum(abs(err{haha})-abs(err{hihi})>0))/64;
    end
end

%
hi=2
err{5}=r_per(:,:,hi)'-p_range(hi);
err{6}=r_BMS(:,:,hi)'-p_range(hi);
err{3}=r_PMP(:,:,hi)'-p_range(hi);
err{4}=r_PMP2(:,:,hi)'-p_range(hi);
err{1}=r_ME(:,:,hi)'-p_range(hi);
err{2}=r_ME2(:,:,hi)'-p_range(hi);
errM2=zeros(6,6);
for haha=1:6
    for hihi=1:6
        errM2(haha,hihi)=sum(sum(abs(err{haha})-abs(err{hihi})>0))/64;
    end
end


