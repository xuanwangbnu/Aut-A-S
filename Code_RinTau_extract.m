clear all;clc

Fidx{1}=dir('Zhenfeng*.mat');
%%
for kk=1:length(Fidx)
    for k=1:length(Fidx{kk})
        %% loading data
        filename=Fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(Fidx{kk})),'--',num2str(kk),'/',num2str(length(Fidx)),'--',filename])
        %% ------------- reading data ---------------------------
        Vm=Ch3.values;dtV=Ch3.interval;     
        Im=Ch4.values;dtI=Ch4.interval;
        Im=Im-mean(Im(1:10));
        if Ch3.units=='pA'
            Im=Im*20;
            Vm=Vm/20;
        end
        
        %%  finding key time point
        TimeStim = eventedgefinderr(Im,-2,1/dtI,0.2,0.5,-1,0);
        TimeStim=TimeStim(:,1).*dtI;

        %% waveform extract
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
            
            Vrest=[Vrest;mean(Vm(baseidxV))];
            waveformV=[waveformV,Vm(dataidxV)];
            waveformI=[waveformI,Im(dataidxI)-mean(Im(baseidxI))];
        end
        tspanV=linspace(-pretime,postime,size(waveformV,1));
        tspanI=linspace(-pretime,postime,size(waveformI,1));
        
        % exclude wrong trace(in case of epileptic network activity)
        wrongidx=[];
        temp=mean(waveformV(tspanV<0|tspanV>0.5,:),1);
        wrongidx=[wrongidx;find(temp>-45)'];
        
        temp=range(waveformV(tspanV>0.6,:),1);
        wrongidx=[wrongidx;find(temp>6)'];
        
        temp=abs(mean(waveformV(tspanV<0,:),1)-mean(waveformV(tspanV>0.6,:),1));
        wrongidx=[wrongidx;find(temp>20)']
        
        temp=range(waveformV(tspanV>0.1&tspanV<0.4,:),1);
        wrongidx=[wrongidx;find(temp>8)'];
        
        temp=min(waveformV(tspanV<0|tspanV>0.5,:),[],1);
        wrongidx=[wrongidx;find(temp<-100)'];
        
        waveformI(:,wrongidx)=[];
        waveformV(:,wrongidx)=[];
        Vrest(wrongidx)=[];
        
        
        %% make same necessary Tag
        % spk shape
        figure(1),clf
        Rin=[];
        tau=[];
        sagRatio=[];
        
        Fun=@(tau,A,B,x) A.*exp(-x./tau)+B;
        fo= fitoptions('Method','NonlinearLeastSquares',...
            'Lower',[0,0.1,-100],...
            'Upper',[10,100,10], ...
            'Startpoint',[0.1,15,-70]);
        
        for i=1:size(waveformV,2)
            x=tspanV(tspanV>0&tspanV<0.1);
            y=waveformV(tspanV>0&tspanV<0.1,i);
            model=fit(x(:),y(:),Fun,fo);
            yy=feval(model,x);
            tau=[tau;(model.tau).*1000];
            
            Max=Vrest(i)-min(waveformV(tspanV>0&tspanV<0.2,i));
            SS=Vrest(i)-mean(waveformV(tspanV>0.4&tspanV<0.45,i));
            sagRatio=[sagRatio;Max./SS];
            
            SSI=-mean(waveformI(tspanI>0.4&tspanI<0.45,i));
            Rin=[Rin;SS./SSI*1000];
            
            figure(1)
            plot(tspanV,waveformV(:,i)),hold on
            plot(x,yy,'r')
            drawnow
            
        end
        
        % StimAmp
        StimAmp=mean(waveformI(tspanI>0.1&tspanI<0.4,:))-mean(waveformI(tspanI<0,:));
        StimAmp=round(StimAmp/10)*10;
 
        save(['DATA_RinTau_',filename],'tspanV','tspanI','waveformV','waveformI','TimeStim',...
            'StimAmp','Vrest','Rin','tau','sagRatio')
    end
    
end

