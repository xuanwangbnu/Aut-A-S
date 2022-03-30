clear;clc;

fidx{1}=dir('Results*Network_isolated0_DA1_dMSNAR0_dMSNAutapase0_iMSNAR0_iMSNAutapase0_FSAR0_FSAutapase0_trial#*.mat');
fidx{2}=dir('Results*Network_isolated0_DA1_dMSNAR0_dMSNAutapase0_iMSNAR0_iMSNAutapase0_FSAR0_FSAutapase1_trial#*.mat');
fidx{3}=dir('Results*Network_isolated0_DA1_dMSNAR0_dMSNAutapase0_iMSNAR0_iMSNAutapase0_FSAR1_FSAutapase1_trial#*.mat');
groupName={'No aut','Aut,SR','Aut,SR+AR'};



%%
for k=1:length(fidx)
    for kk=1:length(fidx{k})
        
        filename=fidx{k}(kk).name;
        load(filename)
        dt=tspan(2)-tspan(1);
        disp([num2str(k),'/',num2str(length(fidx)),'--',filename])
        %%
        for jj=1
            if jj==1
                data=LFPRecord;
            elseif jj==2
                data=LFP_disect_Record(:,1);% FS AR
            elseif jj==3
                data=LFP_disect_Record(:,2);% FS SR
            elseif jj==4
                data=LFP_disect_Record(:,3);% dMSN SR   
            elseif jj==5
                data=LFP_disect_Record(:,4); % iMSN SR                 
            end
            fs=400;
            data=downsample(data,round(1000/dt/fs));
            data=zscore(data);
            [y,f] = PowerSpectrum(data,fs,[0,fs/2],0,0,0);y=y';
            y=y./trapz(f,y);
            y=10*log10(y);
            paramRecord{k}{jj}(:,kk)=y;
        end
        
    end
end
paramName={'LFP','AR','SR','dMSNR','iMSNR'};
save('Results_PSD.mat','groupName','paramName','paramRecord','f')
%%
clear all;
load('Results_PSD.mat')
colorma=lines;
figure(99),clf
for jj=1    
    for k=1:length(groupName)
        data=paramRecord{k}{jj};
        ymean=mean(data,2);
        ysem=std(data,[],2)./sqrt(size(data,2));
        
        plot(f,ymean,'color',colorma(k,:)),hold on
    end
    
    title(paramName{jj})
    ylabel('Power spectral density')
    xlim([0 90])
end
