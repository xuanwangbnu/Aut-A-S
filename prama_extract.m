clear;clc;

fidx=dir('Results*.mat');
fs=100000;
theta=[6 8];
alfa=[8 13];
beta=[13 30];
lowgama=[30 75];
highgama=[75 85];
spikefreq=nan(1,3);
Burstrate=nan(1,3);
Burstnumber=nan(1,3);
Timeinterspike=nan(1,3);
Burstdura=nan(1,3);
Burstpercen=nan(1,3);
Intertrain=nan(1,3);



for k=1:length(fidx)
    filename=fidx(k).name;
    load(filename)
    disp([num2str(k),'/',num2str(length(fidx)),'--',filename])
    %%
    for i=1:length(Ntypes)
        if ~isempty(raster{i})
            spiketime=raster{i}(:,1);
            spiketime=spiketime(spiketime>=1000);
            spikefreq(i)=length(spiketime)./Ntypes(i)./3;
            
            burstrate=nan(Ntypes(i),1);
            burstnumber=nan(Ntypes(i),1);
            Tinterspike=nan(Ntypes(i),1);
            burstdura=nan(Ntypes(i),1);
            intertrain=nan(Ntypes(i),1);
            for j=1:Ntypes(i)
                idx=find(raster{i}(:,2)==j&raster{i}(:,1)>=1000);
                if ~isempty(idx)
                    burstime=raster{i}(idx,1);
                    idx=find(diff(burstime)<16);
                else burstrate(j,1)=0;
                    burstnumber(j,1)=0;
                    Tinterspike(j,1)=0;
                    burstdura(j,1)=0;
                end
                if ~isempty(idx)
                    ind=find(diff(idx)>1);
                    if ~isempty(ind)
                        burstimeR=burstime([idx(ind)+1;idx(end)+1]);
                        burstimeL=burstime([idx(1);idx(ind+1)]);
                        spikenumber=[idx(ind)+1;idx(end)+1]-[idx(1);idx(ind+1)]+2;
                        
                    else burstimeR=burstime(idx(end)+1);
                        burstimeL=burstime(idx(1));
                        spikenumber=[idx(end)+1]-[idx(1)]+2;
                    end
                    interspike=(burstimeR-burstimeL)./(spikenumber-1);
                    burstrate(j,1)=length(burstimeR)./3;%Hz
                    burstnumber(j,1)=mean(spikenumber);
                    Tinterspike(j,1)=mean(interspike);%ms
                    burstdura(j,1)=mean(burstimeR-burstimeL);%ms
                else burstrate(j,1)=0;
                    burstnumber(j,1)=0;
                    Tinterspike(j,1)=0;
                    burstdura(j,1)=0;
                end
                if burstdura(j,1)~=0&length(burstimeL)>1
                    intertrain(j,1)=mean(burstimeL(2:end)-burstimeR(1:end-1));
                else intertrain(j,1)=0;
                end
            end
            
            Burstrate(i)=mean(burstrate);
            Burstnumber(i)=mean(burstnumber);
            Timeinterspike(i)=mean(Tinterspike);
            if max(burstdura)==0
                Burstdura(i)=0;
            else
                Burstdura(i)=mean(burstdura(burstdura~=0));
            end
            ine=find(burstdura~=0);
            Burstpercen(i)=length(ine)./Ntypes(i);
            if max(intertrain)==0
                Intertrain(i)=0;
            else
                Intertrain(i)=mean(intertrain(intertrain~=0));
            end
        else
            Burstrate(i)=0;
            Burstnumber(i)=0;
            Timeinterspike(i)=0;
            Burstdura(i)=0;
            spikefreq(i)=0;
            Burstpercen(i)=0;
            Intertrain(i)=0;
        end
    end
    data=LFPRecord;
    
    fs=400;
    dt=tspan(2)-tspan(1);
    data=downsample(data,round(1000/dt/fs));
    data=zscore(data);
    [y,fspan] = PowerSpectrum(data,fs,[0,fs/2],0,0,0);y=y';
    
    
    % normalized to be power spectrum density
    y=y./trapz(fspan,y);
    
    inx=find(fspan<=alfa(2)&fspan>alfa(1));
    alfapower=trapz(fspan(inx),y(inx));
    inx=find(fspan<=beta(2)&fspan>beta(1));
    betapower=trapz(fspan(inx),y(inx));
    inx=find(fspan<=lowgama(2)&fspan>lowgama(1));
    lowgamapower=trapz(fspan(inx),y(inx));
    inx=find(fspan<=highgama(2)&fspan>highgama(1));
    highgamapower=trapz(fspan(inx),y(inx));
    inx=find(fspan<=theta(2)&fspan>theta(1));
    thetapower=trapz(fspan(inx),y(inx));
    
    save(['Prama_extract_',filename],'spikefreq','Burstrate','Intertrain','Burstnumber','Timeinterspike','Burstdura','Burstpercen','alfapower','betapower','lowgamapower','highgamapower','thetapower')
    
end