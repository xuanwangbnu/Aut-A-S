clear;clc;

fidx=dir('ResultsGD4FS01_Network_isolated0_DA1_dMSNAR0_dMSNAutapase0_iMSNAR0_iMSNAutapase0_FSAR1_FSAutapase1_trial#*.mat');

colorma=lines;
colorma(1,:)=[];

figure(100),clf
wtsum=nan(512,801,length(fidx));
ysum=[];

ccRecordPV2iMSN=[];
ccRecordPV2dMSN=[];
for k=1:length(fidx)
    
    filename=fidx(k).name;
    load(filename)
    disp([num2str(k),'/',num2str(length(fidx)),'--',filename])
    %%
    figure(100+k),clf
    figure(200+k),clf
    figure(300+k),clf
    figure(400+k),clf
    
    ifrRecord=[];
    for i=1:length(Ntypes)
        figure(100+k)
        h1(i)=subplot(length(Ntypes),1,i);
        for j=1:3:size(vRecord{i},2)
            plot(tspan,vRecord{i}(:,j)+(j-1)*0,'color',colorma(i,:)),hold on
        end
        axis tight
        title(filename,'interpreter','none')
        xlim([1000 4000])
        ylabel('Membrane potential')
        figure(200+k)
        h2(i)=subplot(length(Ntypes),1,i);
        if ~isempty(raster{i})
            plot(raster{i}(:,1),raster{i}(:,2),'.','markeredgecolor',colorma(i,:),'markersize',5)
        end
        ylim([0,Ntypes(i)])
        title(filename,'interpreter','none')
        xlim([1000 4000])
        ylabel('AP raster')
        figure(300+k);
        fsFR=10;dt=tspan(2)-tspan(1);
        tt=downsample(tspan,round(1000/dt/fsFR))';
        
        if ~isempty(raster{i})
            ifr=[];
            for j=1:Ntypes(i)
                idx=find(raster{i}(:,2)==j);
                if ~isempty(idx)
                    t_spk=raster{i}(idx,1)/1000;
                    ifr(:,j)=hist(t_spk,tt./1000)/((tt(2)-tt(1))/1000);
                else
                    ifr(:,j)= zeros(size(tt));
                end
            end
            ifrRecord=[ifrRecord,mean(ifr,2)];
            
            Sem=std(ifr,0,2)./sqrt(Ntypes(i));
            Mean=mean(ifr,2);
            ttfill=[tt(:);flipud(tt(:))];
            yfill=[Mean(:)+Sem(:);flipud(Mean(:)-Sem(:))];
            
            h3(i)=subplot(length(Ntypes),1,i);
            fill(ttfill,yfill,colorma(i,:),'edgecolor','none');hold on
            plot(tt,Mean,'color',colorma(i,:))
            axis tight
            
        else
            h3(i)=subplot(length(Ntypes),1,i);
        end
        title(filename,'interpreter','none')
        xlim([1000 4000])
        ylabel('Spike frequency')
        
        figure(400+k)
        h4(i)=subplot(length(Ntypes),1,i);
        data=vMeanRecord(:,i);
        fs=400;
        data=downsample(data,round(1000/dt/fs));
        tt=downsample(tspan,round(1000/dt/fs));
        
        fb=2;fc = 2;
        fspan = linspace(0.005,0.5,512)* fs;
        Scales = fc./ (fspan./fs);
        wt=abs(cwt(data,Scales,['cmor',num2str(fb),'-',num2str(fc)]));
        
        imagesc(tt,fspan,wt);
        caxis([0.01,20])
        set(gca,'ydir','normal')
        
        yyaxis right
        plot(tt,data,'w')
        xlim([1000,4000])
        title(filename,'interpreter','none')
        ylabel('freq_Vpower')
    end
    xlim([0 100])
    linkaxes([h2,h3,h4,h1],'x')
    
    
    figure(500+k)
    data=vMeanRecord(:,3);
    fs=2000;
    data=resample(data,1,round(1000/dt/fs));
    tt=downsample(tspan,round(1000/dt/fs));
    
    fb=4;fc = 2;
    fspan = linspace(0.005,0.5,512)* fs;
    Scales = fc./ (fspan./fs);
    wt=abs(cwt(data,Scales,['cmor',num2str(fb),'-',num2str(fc)]));
    imagesc(tt,fspan,wt)
    caxis([0.01,150])
    set(gca,'ydir','normal')
    ylim([0,100])
    yyaxis right
    plot(tt,data,'w')
    xlim([1000,4000])
    title(filename,'interpreter','none')
    ylabel('LFP time power&LFP trace')
    
    figure(100)
    
    data=vMeanRecord(:,3);
    fs=400;
    data=downsample(data,round(1000/dt/fs));
    
    [y,f] = PowerSpectrum(data,fs,[0,fs/2],0,0,0);y=y';
    y=y./trapz(f,y);
    y=10*log10(y);
    
    plot(f,y),hold on
    ysum=[ysum;y];
    ylabel('Power')
    xlim([0 100])
    
end
%%
figure(99),clf
ymean=mean(ysum,1);
plot(f,ymean),hold on
ylabel('sum_Power')
xlim([0 100])