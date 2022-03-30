clear,clc
fid{1}=dir('*NSr*20060.mat');
fid{2}=dir('*MSr*20060.mat');
for kk=1:length(fid)
preeventfreq=[];
poseventfreq=[];
for k=1:length(fid{kk})
    filename=fid{kk}(k).name;
    load(filename,'prefreq','posfreq');
    preeventfreq=[preeventfreq;nanmean(prefreq)];
    poseventfreq=[poseventfreq;nanmean(posfreq)];    
end
freq=poseventfreq-preeventfreq;
x=kk;
X=kk+0.3;
y=mean(freq);
n=length(freq);
sem=std(freq,[],1)/sqrt(n);
errorbar(x,y,sem,'o','color','k','markerfacecolor','k','markersize',10,'CapSize',2,'LineWidth',1),hold on
plot(X,freq,'ok','markersize',12)
ylabel('deltaEvent freq')
ylim([-10,10])
xlim([0,3])
% A=preeventfreq;
% B=poseventfreq;
% if swtest(A,0.05)==0&&swtest(B,0.05)==0
%     [~,p]=ttest(A,B);
%     disp('ttest')
% else
%     p=signrank(A,B);
%     disp('signrank')
% end
% text(1,10,num2str(p))
end

