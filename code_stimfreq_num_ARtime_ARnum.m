clear,clc
fid=dir('*Result_200*.mat');
totalARtime=[];
totalARnum=[];
totalfreq=[];
totalnum=[];
colorma=lines;
for k=1:length(fid)
    filename=fid(k).name;
    load(filename);
    totalARtime=[totalARtime,ARtime];
    totalARnum=[totalARnum,ARnum];
    totalfreq=[totalfreq,freq];
    totalnum=[totalnum,num];
    totalnum=round(totalnum/20)*20;
end
uninum=unique(totalnum);
uninum=uninum(uninum<70);
numARtime=nan(length(fid),length(uninum));
numARnum=nan(length(fid),length(uninum));
for i=1:length(uninum)    
    idx=find(totalnum==uninum(i));
    numARtime(1:length(idx),i)=totalARtime(idx);
    numARnum(1:length(idx),i)=totalARnum(idx);
end

figure(1),clf
x=uninum;
yt=nanmean(numARtime,1);
nt=sum(~isnan(numARtime),1);
st=nanstd(numARtime,[],1)./sqrt(nt);
errorbar(x,yt,st,'-o','color','k','markerfacecolor','w','markersize',6,'CapSize',6,'LineWidth',1),hold on
plot(x,numARtime,':k','LineWidth',0.5)
ylabel('PT-AR duration(ms)')
xlabel('APs')
xlim([10 70])
title('numARtime')

figure(2),clf
yn=nanmean(numARnum,1);
nn=sum(~isnan(numARnum),1);
sn=nanstd(numARnum,[],1)./sqrt(nn);
errorbar(x,yn,sn,'-o','color','k','markerfacecolor','w','markersize',6,'CapSize',6,'LineWidth',1),hold on
plot(x,numARnum,':k','LineWidth',1)
  ylabel('PT-AR events')
xlabel('APs')
xlim([10 70])
title('numARnum')

PTARtime=yt
PTARtimest=st
PTARnum=yn
PTARnumst=sn
NUM=numARnum(~isnan(numARnum));
NUM=reshape(NUM,[],3);
TIME=numARtime(~isnan(numARtime));
TIME=reshape(TIME,[],3);

groupNUM=[ones(size(NUM,1),1);2*ones(size(NUM,1),1);3*ones(size(NUM,1),1)];
[pNUM,tNUM,statsNUM]=kruskalwallis(NUM(:),groupNUM,'off');
cNUM = multcompare(statsNUM,'Display','off');
groupTIME=[ones(size(TIME,1),1);2*ones(size(TIME,1),1);3*ones(size(TIME,1),1)];
[pTIME,tTIME,statsTIME]=kruskalwallis(TIME(:),groupTIME,'off');
cTIME = multcompare(statsTIME,'Display','off');
