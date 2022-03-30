clear,clc
fid{1}=dir('*NS*.mat');
fid{2}=dir('*Sr*.mat');
for kk=1:length(fid)
    for k=1:length(fid{kk})
        filename=fid{kk}(k).name;
        load(filename);
        if kk==1
            ctrlARtime(k)=ARtime;
            ctrlARnum(k)=ARnum;
        else
            SrARtime(k)=ARtime;
            SrARnum(k)=ARnum;
        end
    end
end
clen=length(ctrlARtime);
slen=length(SrARtime);
if clen>=slen
    SrARtime(slen+1:clen)=nan;
    SrARnum(slen+1:clen)=nan;
else
    ctrlARtime(clen+1:slen)=nan;
    ctrlARnum(clen+1:slen)=nan;
end
x=1:2;
yt=[ctrlARtime',SrARtime'];
yn=[ctrlARnum',SrARnum'];
Mt=nanmean(yt,1);
Mn=nanmean(yn,1);
Nt=sum(~isnan(yt),1);
Nn=sum(~isnan(yn),1);
St=nanstd(yt,[],1)./sqrt(Nt);
Sn=nanstd(yn,[],1)./sqrt(Nn);
figure(1),clf
bar(x,Mt,0.14,'FaceColor',[1,1,1],'EdgeColor','k','lineWidth',1),hold on
errorbar(x,Mt,St,'o','color','k','markersize',1,'CapSize',3,'LineWidth',1)
plot(x-0.05+0.1*rand(size(yt,1),1),yt,'ob','MarkerSize',10.5)
ylabel('PT-AR duration(ms)')
xlim([0.5,2.5])
A=ctrlARtime(~isnan(ctrlARtime));
B=SrARtime(~isnan(SrARtime));
if swtest(A,0.05)==0&&swtest(B,0.05)==0
    [~,p]=ttest2(A,B);
    display('ttest')
else
    p=ranksum(A,B);
    display('ranksum')
end
text(1,100,num2str(p))
figure(2),clf
bar(x,Mn,0.14,'FaceColor',[1,1,1],'EdgeColor','k','lineWidth',1),hold on
errorbar(x,Mn,Sn,'o','color','k','markersize',1,'CapSize',3,'LineWidth',1)
plot(x-0.05+0.1*rand(size(yn,1),1),yn,'ob','MarkerSize',10.5)
ylabel('PT-AR number')
xlim([0.5,2.5])
A=ctrlARnum(~isnan(ctrlARnum));
B=SrARnum(~isnan(SrARnum));
if swtest(A,0.05)==0&&swtest(B,0.05)==0
      display('ttest')
    [~,p]=ttest2(A,B);
else
    p=ranksum(A,B);
     display('ranksum')
end
text(1,10,num2str(p))
PTARtime=Mt
PTARtimest=St
PTARnum=Mn
PTARnumst=Sn

