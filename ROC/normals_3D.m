%% 3D / mean shift

% colorM0=colormap('autumn'); % Get color matrix
% color_id=floor(linspace(1,256,length(Ns_range)));
% colorM=colorM0(color_id,:);
Ns_range=10:10:200;
mu_range=0.5:0.25:3;
p_range=[.3 .7];
r_BMS=zeros(length(Ns_range),length(mu_range),length(p_range));
r_ME=zeros(length(Ns_range),length(mu_range),length(p_range));
r_ME2=zeros(length(Ns_range),length(mu_range),length(p_range));
r_PMP=zeros(length(Ns_range),length(mu_range),length(p_range));
r_PMP2=zeros(length(Ns_range),length(mu_range),length(p_range));
r_per=zeros(length(Ns_range),length(mu_range),length(p_range));
for haha=1:length(p_range)
    haha
    for s=1:length(Ns_range)
        fprintf('|||')
        for mu=1:length(mu_range)
            % sample normal distributions
            Ndata=100;
            Nsample=Ns_range(s); % N trial
            p_1=p_range(haha);
            data1=zeros(Nsample, ceil(Ndata*p_1));
            data2=zeros(Nsample, Ndata-ceil(Ndata*p_1));
            for i=1:ceil(Ndata*p_1)
                data1(:,i)=2*randn(Nsample,1)+mu_range(mu);
            end
            for i=1:Ndata-ceil(Ndata*p_1)
                data2(:,i)=2*randn(Nsample,1);
            end
            data_all=[data1, data2];
            % fit normal distributions
            norm2pdf=@(data,sd)normpdf(data,0,sd);
            param1=zeros(Ndata,2);
            for i=1:Ndata
                param1(i,:)=mle(data_all(:,i),'distribution','normal', 'lowerbound', [-10. 0.001], 'upperbound', [10, 100]);
            end
            param2=zeros(Ndata,1);
            for i=1:Ndata
                param2(i,:)=mle(data_all(:,i),'pdf',norm2pdf, 'start', 0.5, 'lowerbound', 0.001, 'upperbound', 100);
            end
            % get BICs
            BIC_1=zeros(Ndata,1);
            BIC_2=zeros(Ndata,1);
            for i=1:Ndata
                BIC1=log(Nsample)*2;
                BIC2=log(Nsample)*1;
%                 BIC1=4;
%                 BIC2=2;
                for j=1:Nsample
                    BIC1=BIC1-2*log(normpdf(data_all(j,i),param1(i,1),param1(i,2)));
                    BIC2=BIC2-2*log(normpdf(data_all(j,i),0,param2(i)));
                end
                BIC_1(i)=BIC1;
                BIC_2(i)=BIC2;
            end
            BIC_all0=[BIC_1,BIC_2]';
%             
%             plot(BIC_all(2,:)-BIC_all(1,:),'Color',colorM(s, :))
%             hold on
%             223,124,1005,651
            BIC_all=BIC_all0-repmat(min(BIC_all0,[],1),[2,1]);
            % RFX-BMS
            Opt_BMC.Start=1e-6*ones(2,1); % Flat prior
            Opt_BMC.MaxIter=1e6;
            Opt_BMC.Stop=1e-6;
            Opt_BMC.Verbosity=0;
            Opt_BMC.Rec=0;
            BMS_all=BMW_BMS(-BIC_all, Opt_BMC);
            r_BMS(s,mu,haha)=BMS_all.r(1);
            % mean/median model evidence
            BIC_mean=mean(BIC_all0,2);
            BIC_mean=BIC_mean-repmat(min(BIC_mean,[],1),[2,1]);
            BIC_mean2=log(mean(exp(BIC_all),2));
            r_ME(s,mu,haha)=exp(-.5*BIC_mean(1,:))./sum(exp(-.5*BIC_mean));
            r_ME2(s,mu,haha)=exp(-.5*BIC_mean2(1,:))./sum(exp(-.5*BIC_mean2));
            % mean/median posterior model probability, p(model | data)
            pmp_1=exp(-.5*BIC_all(1,:))./sum(exp(-.5*BIC_all));
            pmp_2=exp(-.5*BIC_all(2,:))./sum(exp(-.5*BIC_all));
            r_PMP(s,mu,haha)=mean(pmp_1);
            r_PMP2(s,mu,haha)=median(pmp_1);
            % percentage of being the best model
            r_per(s,mu,haha)=mean(BIC_1<=BIC_2);
        end
    end
end

% Visualization
% colorM = [0.01*ones(length(Ns_range),1),(0.1:(0.9/length(Ns_range)):(1-0.9/length(Ns_range)))',  0.01*ones(length(Ns_range),1)]; 
colorM0=colormap('autumn'); % Get color matrix
color_id=floor(linspace(1,256,length(Ns_range)));   
colorM=colorM0(color_id,:);
hi=2;
fighaha=figure(1);
set(fighaha,'Position',[71,73,1500,756])
for hei=1:length(Ns_range)
    %---%
    subplot(3,2,6)
    title('RFX-BMS')
    if all(abs(r_BMS(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 12])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticks('')
    plot(r_BMS(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(mu_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,2)
    title('Mean PMP')
    if all(abs(r_PMP(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 12])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticks('')
    plot(r_PMP(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(mu_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,4)
    title('Median PMP')
    if all(abs(r_PMP2(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 12])
    xticks('')
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticklabels({num2str(mu_range(1)), num2str(mu_range(round(length(mu_range)/2))), num2str(mu_range(end))})
    plot(r_PMP2(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(mu_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,5)
    title('Best Model Percentage')
    xlim([0 12])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xlabel('Mu Shift')
    xticks([1, round(length(mu_range)/2), length(mu_range)])
    if all(abs(r_per(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    plot(r_per(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(mu_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,1)
    title('Mean Log Model Evidence')
    ylabel('Error of Recovery')
    if all(abs(r_ME(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 12])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticks('')
    plot(r_ME(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(mu_range)],[0 0],'Color',[0,0,0])
    %---%
    subplot(3,2,3)
    title('Mean Model Evidence')
    if all(abs(r_ME2(hei,:,hi)'-p_range(hi))<=0.3), ylim([-0.3,0.3]), end
    xlim([0 12])
    ylim([0-p_range(hi), 1-p_range(hi)])
    xticks('')
    plot(r_ME2(hei,:,hi)'-p_range(hi), '.-','Color',colorM(hei, :))
    hold on
    plot([0, 1+length(mu_range)],[0 0],'Color',[0,0,0])
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
        errM1(haha,hihi)=sum(sum(abs(err{haha})-abs(err{hihi})>0))/220;
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
        errM2(haha,hihi)=sum(sum(abs(err{haha})-abs(err{hihi})>0))/220;
    end
end


