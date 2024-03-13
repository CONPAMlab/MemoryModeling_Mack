%% Prediction Check
subj_pred_id = find(1-FitResults_new{2}.Param(:,end)<= 0.0001);
param_DPSD_fit = FitResults_new{1}.Param(subj_pred_id,:);
param_DPSD = [param_DPSD_fit(:,1), -100*ones(length(subj_pred_id),1), param_DPSD_fit(:,2:(end-1)), 100*ones(length(subj_pred_id),1), ones(length(subj_pred_id),1), param_DPSD_fit(:,end), zeros(length(subj_pred_id),1)];
%[d-prime, c, sigma, Ro, Rn]
param_Slot_fit = FitResults_new{2}.Param(subj_pred_id,:);
param_Slot = [param_Slot_fit(:,1), -100*ones(length(subj_pred_id),1), param_Slot_fit(:,2:(end-1)), 100*ones(length(subj_pred_id),1), ones(length(subj_pred_id),1), param_Slot_fit(:,end), zeros(length(subj_pred_id),1)];
%[d-prime, c, sigma, Pm, Pn]
Phit_pred = zeros(length(subj_pred_id),7);
Pfa_pred = zeros(length(subj_pred_id),7);
Phit_pred2 = zeros(length(subj_pred_id),7);
Pfa_pred2 = zeros(length(subj_pred_id),7);
for i = 1:length(subj_pred_id)
%     [data_sim, ~]=SimDPSD(param_DPSD(i,:), 5000, 5000, 1, 5);
    [data_sim, ~]=SimMixture(param_Slot(i,:), 5000, 5000, 1, 7);
    Phit_pred(i,:) = data_sim{1}.Nhit(:,2:end)/5000;
    Pfa_pred(i,:) = data_sim{1}.Nfa(:,2:end)/5000;
    Phit_real(i,:) = DataHAHAMerged{subj_pred_id(i)}.Nhit(2:end)/DataHAHAMerged{subj_pred_id(i)}.Nhit(1);
    Pfa_real(i,:) = DataHAHAMerged{subj_pred_id(i)}.Nfa(2:end)/DataHAHAMerged{subj_pred_id(i)}.Nfa(1);
end
for i = 1:length(subj_pred_id)
    [data_sim, ~]=SimDPSD(param_DPSD(i,:), 5000, 5000, 1, 7);
%     [data_sim, ~]=SimMixture(param_Slot(i,:), 5000, 5000, 1, 5);
    Phit_pred2(i,:) = data_sim{1}.Nhit(:,2:end)/5000;
    Pfa_pred2(i,:) = data_sim{1}.Nfa(:,2:end)/5000;
    Phit_real(i,:) = DataHAHAMerged{subj_pred_id(i)}.Nhit(2:end)/DataHAHAMerged{subj_pred_id(i)}.Nhit(1);
    Pfa_real(i,:) = DataHAHAMerged{subj_pred_id(i)}.Nfa(2:end)/DataHAHAMerged{subj_pred_id(i)}.Nfa(1);
end

for i = 1:length(subj_pred_id)
    plot(Pfa_pred(i,:), Phit_pred(i,:),'-')
    hold on
    plot(Pfa_pred2(i,:), Phit_pred2(i,:),'-')
    plot(Pfa_real(i,:), Phit_real(i,:),'o')
    xlim([0 1])
    ylim([0 1])
    xlabel('P(fa)')
    ylabel('P(hit)')
    pause
    hold off
end