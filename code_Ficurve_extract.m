clear all;clc
Fidx{1}=dir('Zhenfeng*.mat');
for kk=1:length(Fidx)
        for k=1:length(Fidx{kk})
        
        %loading data
        filename=Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename]);
        %%
        %reading data
        Vm=Ch3.values;dtV=Ch3.interval;     
        Im=Ch4.values;dtI=Ch4.interval;
         Im=Im-mean(Im(1:10));
        if Ch3.units=='pA'
            Im=Im*20;
            Vm=Vm/20;
        end
        
        %finding key time point
        TimeStim = eventedgefinderr(Im,100,1/dtI,0.2,0.5,1,0).*dtI;
        TimeStim=TimeStim(:,1);
        disp(length(TimeStim));
        
        %waveform extract
        pretime=0.1;
        postime=0.8;
        
        waveformV=[];
        waveformI=[];
        Vrest=[];
        for i=1:length(TimeStim)
            dataidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV+postime/dtV);
            baseidxV=round(TimeStim(i)/dtV-pretime/dtV:TimeStim(i)/dtV);
            dataidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI+postime/dtI);
            baseidxI=round(TimeStim(i)/dtI-pretime/dtI:TimeStim(i)/dtI);
            
            Vrest=[Vrest,mean(Vm(baseidxV))];
            waveformV=[waveformV,Vm(dataidxV)];
            waveformI=[waveformI,Im(dataidxI)-mean(Im(baseidxI))];
        end
        tspanV=linspace(-pretime,postime,size(waveformV,1));
        tspanI=linspace(-pretime,postime,size(waveformI,1));
        
        wrongidx=[];
        temp=mean(waveformV(tspanV<0|(tspanV>0.6&tspanV<1),:),1);
        wrongidx=[wrongidx;find(temp>-45)'];
        temp=mean(waveformV((tspanV>0.6),:),1);
        wrongidx=[wrongidx;find(temp>-50)'];
       temp=range(waveformV((tspanV>0.6),:));
       wrongidx=[wrongidx;find(temp>12)'];
        temp=min(waveformV(tspanV<0|tspanV>0.5,:),[],1);
        wrongidx=[wrongidx;find(temp<-100)'];
        
        temp=mean(waveformV(tspanV<0,:),1);
        wrongidx=[wrongidx;find(temp>-55)']
        
        
        waveformI(:,wrongidx)=[];
        waveformV(:,wrongidx)=[];
        Vrest(wrongidx)=[];
        TimeStim(wrongidx)=[];
        
        %make same necessary Tag
        %spk shape
        spkFreq=[];
        spkTime=[];
        
        for i=1:size(waveformV,2)
            peaktime=peakfinder(waveformV(:,i)-Vrest(i),10,0, 1, 0).*dtV;
            worngidx=[];
            for j=1:length(peaktime)
                idx=round(peaktime(j)/dtV-0.01/dtV:peaktime(j)/dtV);
                temp=min(waveformV(idx,i));
                if temp>waveformV(idx(end),i)-10
                    worngidx=[worngidx;j];
                end
                
                if waveformV(idx(end),i)<-20&peaktime(j)>0.5+pretime
                    worngidx=[worngidx;j];
                end
                
            end
            peaktime(worngidx)=[];
            
            if ~isempty(peaktime)
                peaktime(peaktime>0.51+pretime|peaktime<pretime+0.0005)=[];
            end
        
            spkFreq=[spkFreq;length(peaktime)/0.5];
            spkTime{i}=peaktime;
            
            figure(1),clf
            plot(tspanV,waveformV(:,i)),hold on 
            plot(tspanV(round(peaktime/dtV)),waveformV(round(peaktime/dtV),i),'ro')
            drawnow
        end
        
        % StimAmp
        StimAmp=mean(waveformI(tspanI>0.1&tspanI<0.4,:))-mean(waveformI(tspanI<0,:));
        
        % PC PV protocol
        StimAmp=round(StimAmp/10)*10;
        
        save(['DATA_FIcurve_',filename],'tspanV','tspanI','waveformV','waveformI','TimeStim',...
            'Vrest','spkFreq','spkTime','StimAmp')
        
    end
    
end