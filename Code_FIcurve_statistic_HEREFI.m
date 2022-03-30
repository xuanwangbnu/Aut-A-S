clear all;clc


Fidx{1}=dir('DATA_FIcurve_*.mat');

GroupName={'MSN'};

Num=2;
thAP=1;

Binwidth=100;
BinL=50:Binwidth:2800-Binwidth;
BinR=BinL+Binwidth;

BinsFIcurve{1}=[BinL',BinR'];





%%

colorma=lines;

VposiRecord=[];

for kk=1:length(Fidx)
    for k=1:length(Fidx{kk})
   
        filename=Fidx{kk}(k).name;
        load(filename)
        dt=tspanV(2)-tspanV(1);
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        FilenameRecord{kk}{k,1}=filename;
        
        %% fitting ReLU
        Fun=@(a,b,x) max([zeros(length(x),1),a.*x(:)+b],[],2);
        fo= fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0.001,-1000],...
            'Upper',[10,1000], ...
            'Startpoint',[0.1,-10]);
        
        uniqueStim=unique(StimAmp);
        StimAmp_act=zeros(size(uniqueStim))';
        spkFreq_act=zeros(size(uniqueStim))';
        for i=1:length(uniqueStim)
            idx=find(StimAmp==uniqueStim(i));
            StimAmp_act(i)=mean(StimAmp(idx));
            spkFreq_act(i)=mean(spkFreq(idx));
        end
        
        [tempx,idx]=sort(StimAmp_act);tempy=spkFreq_act(idx);
        idx=find(tempy>=12,1,'first');
        x=[0;tempx(1:idx)];
        y=[0;tempy(1:idx)];
        
        [Model,~]=fit(x(:),y(:),fittype(Fun),fo);
        xx=linspace(0,x(end),100);
        yy=feval(Model,xx);
      
        FIcurveParamRecord{kk}(k,:)=[-(Model.b/Model.a),Model.a];
        %% ------------- Bin FIcurve ---------------------------
        
        for i=1:size(BinsFIcurve{kk},1)
            idx=find(StimAmp>=BinsFIcurve{kk}(i,1)&StimAmp<BinsFIcurve{kk}(i,2));
            if ~isempty(idx)
                spkFreqBin{kk}(k,i)=mean(spkFreq(idx));
            else
                spkFreqBin{kk}(k,i)=nan;
            end
        end
        
        
        %% ------------- Bin spkshape ---------------------------
        
        Wantidx=find((spkFreq>=Num*2)&(spkFreq>=10*2)&StimAmp(:)<800);
        
        if ~isempty(Wantidx)
            spkwaveform_temp=[];
            HWtemp=[];
            AHPtemp=[];
            Thresholdtemp=[];
            Amptemp=[];
            max_slopetemp=[];
            min_slopetemp=[];
            
            for jj=1:length(Wantidx)
                data=waveformV(:,Wantidx(jj));
                [spkwaveform,ttspanspk,...
                    HW, AHP, Threshold, Amp,max_slope,min_slope,BW,AHPArea,slopeSpecific,...
                    HWidx,AHPidx, Thresholdidx,max_slopeidx,min_slopeidx,BWidx] = AP_Statistic_HERE(data,spkTime{Wantidx(jj)},tspanV,1:max([Num-1,1]));
                
                spkwaveform_temp(jj,:,:)=spkwaveform;
                HWtemp(:,jj)=HW;
                AHPtemp(:,jj)=AHP;
                Thresholdtemp(:,jj)=Threshold;
                Amptemp(:,jj)=Amp;
                max_slopetemp(:,jj)=max_slope;
                min_slopetemp(:,jj)=min_slope;
                
                figure(1),clf
                H=[];
                for j=1:size(spkwaveform,2)
                    H(j)=subplot(1,size(spkwaveform,2),j);
                    plot(ttspanspk,spkwaveform(:,j)),hold on
                    plot(ttspanspk(HWidx(j,:)),spkwaveform(HWidx(j,:),j),'ro')
                    plot(ttspanspk(AHPidx(j)),spkwaveform(AHPidx(j),j),'ko')
                    plot(ttspanspk(Thresholdidx(j)),spkwaveform(Thresholdidx(j),j),'co')
                    axis tight
                    drawnow
                end
                linkaxes(H,'x')
            end
            spkparamRecord{kk}{1}(:,k)=mean(HWtemp,2);
            spkparamRecord{kk}{2}(:,k)=mean(AHPtemp,2);
            spkparamRecord{kk}{3}(:,k)=mean(Thresholdtemp,2);
            spkparamRecord{kk}{4}(:,k)=mean(Amptemp,2);
            spkparamRecord{kk}{5}(:,k)=mean(max_slopetemp,2);
            spkparamRecord{kk}{6}(:,k)=mean(min_slopetemp,2);
            
            if size(spkwaveform_temp,2)==2600
                spkwaveformRecord{kk}(k,:,:)=squeeze(mean(spkwaveform_temp,1));
            else
                spkwaveformRecord{kk}(k,:,:)=nan(2600,Num-1);
            end
        else
            spkparamRecord{kk}{1}(:,k)=nan(Num-1,1);
            spkparamRecord{kk}{2}(:,k)=nan(Num-1,1);
            spkparamRecord{kk}{3}(:,k)=nan(Num-1,1);
            spkparamRecord{kk}{4}(:,k)=nan(Num-1,1);
            spkparamRecord{kk}{5}(:,k)=nan(Num-1,1);
            spkparamRecord{kk}{6}(:,k)=nan(Num-1,1);
            spkwaveformRecord{kk}(k,:,:)=nan(2600,Num-1);
        end
     
    end
end

FIcurveParamName={'theta','slope'};
spkParamName={'HW', 'AHP', 'Threshold', 'Amp','max_slope','min_slope'};

save('Results_FIcurveStatistic.mat','Fidx','GroupName','colorma','tspanV',...
    'spkwaveformRecord','spkparamRecord',...
    'spkFreqBin','BinsFIcurve','FIcurveParamRecord',...
    'Num','thAP',...
    'FIcurveParamName','spkParamName','GroupName','FilenameRecord')

for kk=1:length(Fidx)
    
    xlswrite('Results_FIcurveStatistic.xls',FilenameRecord{kk},['FIcurve_',GroupName{kk}],'A2')
    xlswrite('Results_FIcurveStatistic.xls',FIcurveParamName,['FIcurve_',GroupName{kk}],'B1')
    xlswrite('Results_FIcurveStatistic.xls',FIcurveParamRecord{kk},['FIcurve_',GroupName{kk}],'B2')
    
    xlswrite('Results_FIcurveStatistic.xls',FilenameRecord{kk},['spkshape_',GroupName{kk}],'A2')
    xlswrite('Results_FIcurveStatistic.xls',spkParamName,['spkshape_',GroupName{kk}],'B1')
    xlswrite('Results_FIcurveStatistic.xls',[spkparamRecord{kk}{1}(thAP,:)',...
                                             spkparamRecord{kk}{2}(thAP,:)',...
                                             spkparamRecord{kk}{3}(thAP,:)',...
                                             spkparamRecord{kk}{4}(thAP,:)',...
                                             spkparamRecord{kk}{5}(thAP,:)',...
                                             spkparamRecord{kk}{6}(thAP,:)'],['spkshape_',GroupName{kk}],'B2')


    xlswrite('Results_FIcurveStatistic.xls',FilenameRecord{kk},['FIcurve_',GroupName{kk}],'A2')
    xlswrite('Results_FIcurveStatistic.xls',BinsFIcurve{kk}(:,1)',['FIcurve_',GroupName{kk}],'B1')
    xlswrite('Results_FIcurveStatistic.xls',spkFreqBin{kk},['FIcurve_',GroupName{kk}],'B2')
                  
end

%%

clear all;clc
load('Results_FIcurveStatistic.mat')

% ==================== FIcurve ====================
figure(1),clf
for kk=1:length(Fidx)
    
    x=BinsFIcurve{kk}(:,1);
    y=nanmean(spkFreqBin{kk},1);
    N=sum(~isnan(spkFreqBin{kk}),1);
    sem=nanstd(spkFreqBin{kk},[],1)./sqrt(N);
    
    errorbar(x,y,sem,'r-o','color',colorma(kk,:),'markerfacecolor',colorma(kk,:),'markersize',8),hold on
    
end
legend(GroupName)
ylabel('Frequency (Hz)')
xlabel('Current injection (pA)')
xlim([0 1050])


% ====================  statistic ====================
figure(4),clf
for k=1:length(FIcurveParamName)
    subplot(1,length(FIcurveParamName),k)
    for kk=1:length(Fidx)
        data{kk}=FIcurveParamRecord{kk}(:,k);N=sum(~isnan(data{kk}));
        Mean=nanmean(data{kk});
        Sem=nanstd(data{kk})./sqrt(N);
        bar(kk,Mean,'facecolor',colorma(kk,:));hold on
        errorbar(kk,Mean,Sem,'color',colorma(kk,:))
    end
    ylabel(FIcurveParamName{k})
end

figure(5),clf
for k=1:length(spkParamName)
    
    subplot(1,length(spkParamName),k)
    for kk=1:length(Fidx)
        data{kk}=spkparamRecord{kk}{k}(thAP,:);N=sum(~isnan(data{kk}));
        Mean=nanmean(data{kk});
        Sem=nanstd(data{kk})./sqrt(N);
        bar(kk,Mean,'facecolor',colorma(kk,:));hold on
        errorbar(kk,Mean,Sem,'color',colorma(kk,:))
    end
    ylabel(spkParamName{k}) 
end

% ====================  waveform  ====================

figure(6),clf

xlim_max=0.008;
xlim_min=-0.0015;

for kk=1:length(Fidx)
    for j=1:size(spkwaveformRecord{kk},3)
        subplot(2,size(spkwaveformRecord{kk},3),j);
        data=squeeze(nanmean(spkwaveformRecord{kk}(:,:,j),1));
        plot(ttspanspk,data,'color',colorma(kk,:)),hold on
        axis([xlim_min,xlim_max,-90,60])
        
        subplot(2,size(spkwaveformRecord{kk},3),j+size(spkwaveformRecord{kk},3))
        tempidx=find(ttspanspk>-0.004&ttspanspk<0.0040);
        plot(data(tempidx(2:end)),diff(data(tempidx))*100,'color',colorma(kk,:)),hold on
        axis([-90,60,-1000,1800])
        
    end
end

