clear;clc
Fidx=dir('DATA_RinTau_*.mat');
totalRin=nan(1,length(Fidx));
totalVrest=nan(1,length(Fidx));
for k=1:length(Fidx)
    filename=Fidx(k).name;
    load(filename,'Vrest','Rin','tau','sagRatio','StimAmp','waveformV','tspanV')
    disp([num2str(k),'/','--','/',num2str(length(Fidx)),'--',filename])
    if ~isempty(Rin)&~isempty(Vrest)
    totalRin(k)=mean(Rin);
    totalVrest(k)=mean(Vrest); 
    end
end
MRin=nanmean(totalRin)
NRin=sum(~isnan(totalRin));
SemRin=nanstd(totalRin)/sqrt(NRin)
MVrest=nanmean(totalVrest)
NVrest=sum(~isnan(totalVrest));
SemVrest=nanstd(totalVrest)/sqrt(NVrest)