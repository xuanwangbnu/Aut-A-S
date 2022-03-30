
clear;clc

load('dt.mat');
fid=dir('dataI_*q200_*b60_*.mat');
freq=200;
num=60;
SRtimewindow=0.004;
stimlatency=0.0015  ;
pretime=0.3;
postime=pretime;
%——————————————
bintime=0.05;
arrange=1;
baserange=2;
sti=0:(1/freq):(1/freq)*num;
sti=sti(1:num)';
sti=sti+baserange;
binLbase=[0:bintime:baserange-bintime];
binRbase=[bintime:bintime:baserange];
for k=1:length(fid)

    filename=fid(k).name;
    xlsname=[filename(1:end-4) '.xls'];
    load(filename);
    current=data_act;
    xscale=(1:length(data_act))*dtI;
    [Data,DataText]=xlsread(xlsname);
    total=Data(2:size(Data,1),2);
    total=total./1000;
    binLpt=[sti(end)+stimlatency+SRtimewindow:bintime:round(length(current)*dtI)-bintime];
    binRpt=[sti(end)+stimlatency+SRtimewindow+bintime:bintime:round(length(current)*dtI)];
    binL=[binLbase,binLpt];
    binR=[binRbase,binRpt];
    bin=[binL',binR'];

    for i=1:size(bin,1)
        idx=find(total<bin(i,2)&total>=bin(i,1));
        binARevent(k,i)=length(idx);
        binARfreq(k,i)=binARevent(k,i)/bintime;
    end
    figure(1),clf
    H1=subplot(2,1,1);
    plot(xscale,current)
    xlim([1.4,4.0])
    ylim([-50,20])
    title(filename)
    xlabel('Time(s)')
    ylabel('Current(pA)')
    
    H2=subplot(2,1,2);
    plot(mean(bin,2),binARfreq(k,:),'o-b')
    xlabel('Time(s)')
    ylabel('event frequency(Hz)')
    xlim([1.6,3.2])
    linkaxes([H1 H2],'x')
    drawnow
    pause
    ind=find(total>baserange-pretime&total<baserange);
    eventpre(k)=length(ind);
    ind=find(total>=sti(end)+stimlatency+SRtimewindow&total<=sti(end)+stimlatency+SRtimewindow+postime);
    eventpos(k)=length(ind); 
end

figure(2),clf
x=nanmean(bin,2);
y=nanmean(binARfreq,1);
n=sum(~isnan(binARfreq),1);
s=nanstd(binARfreq,[],1)./sqrt(n);
errorbar(x,y,s,'o','CapSize',1,'MarkerSize',8)
 xlabel('Time(s)')
    ylabel('event frequency(Hz)')
    xlim([1.6,3.2])

time=x;
eventfreq=y;
eventpre=mean(eventpre);
eventpos=mean(eventpos);


                                                                                                                                                                                                                                                                     
