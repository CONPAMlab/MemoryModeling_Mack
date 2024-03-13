Ns_range=10:10:200;
mu_range=0.5:0.25:3;
r_BMS=zeros(length(Ns_range),length(mu_range));
r_ME=zeros(length(Ns_range),length(mu_range));
r_ME2=zeros(length(Ns_range),length(mu_range));
r_PMP=zeros(length(Ns_range),length(mu_range));
r_PMP2=zeros(length(Ns_range),length(mu_range));
r_per=zeros(length(Ns_range),length(mu_range));
for s=1:length(Ns_range)
    for mu=1:length(mu_range)
    % sample normal distributions
    Ndata=100;
    Nsample=Ns_range(s);
    p_1=0.7;
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
        for j=1:Nsample
            BIC1=BIC1-2*log(normpdf(data_all(j,i),param1(i,1),param1(i,2)));
            BIC2=BIC2-2*log(normpdf(data_all(j,i),0,param2(i)));
        end
        BIC_1(i)=BIC1;
        BIC_2(i)=BIC2;
    end
    BIC_all=[BIC_1,BIC_2]';
    BIC_all=BIC_all-repmat(min(BIC_all,[],1),[2,1]);
    % RFX-BMS
    Opt_BMC.Start=1e-6*ones(2,1); % Flat prior
    Opt_BMC.MaxIter=1e6;
    Opt_BMC.Stop=1e-6;
    Opt_BMC.Verbosity=0;
    Opt_BMC.Rec=0;
    BMS_all=BMW_BMS(-BIC_all, Opt_BMC);
    r_BMS(s,mu)=BMS_all.r(1);
    % mean/median model evidence 
    BIC_mean=mean(BIC_all,2);
    BIC_median=median(BIC_all,2);
    r_ME(s,mu)=exp(-.5*BIC_mean(1,:))./sum(exp(-.5*BIC_mean));
    r_ME2(s,mu)=exp(-.5*BIC_median(1,:))./sum(exp(-.5*BIC_median));
    % mean/median posterior model probability, p(model | data)
    pmp_1=exp(-.5*BIC_all(1,:))./sum(exp(-.5*BIC_all));
    pmp_2=exp(-.5*BIC_all(2,:))./sum(exp(-.5*BIC_all));
    r_PMP(s,mu)=mean(pmp_1);
    r_PMP2(s,mu)=median(pmp_1);
    % percentage of being the best model
    r_per(s,mu)=mean(BIC_1<=BIC_2);
    end
end

subplot(3,2,1)
plot(r_BMS'-p_1, 'b.-')
hold on
plot([1, length(mu_range)],[0 0])
subplot(3,2,3)
plot(r_PMP'-p_1, 'b.-')
hold on
plot([1, length(mu_range)],[0 0])
subplot(3,2,5)
plot(r_PMP2'-p_1, 'b.-')
hold on
plot([1, length(mu_range)],[0 0])
subplot(3,2,2)
plot(r_per'-p_1, 'b.-')
hold on
plot([1, length(mu_range)],[0 0])
subplot(3,2,4)
plot(r_ME'-p_1, 'b.-')
hold on
plot([1, length(mu_range)],[0 0])
subplot(3,2,6)
plot(r_ME2'-p_1, 'b.-')
hold on
plot([1, length(mu_range)],[0 0])