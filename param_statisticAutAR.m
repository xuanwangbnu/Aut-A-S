clear;clc;
fidx{1}=dir('Prama_extract_ResultsGD4FS01_Network_isolated0_DA1_dMSNAR0_dMSNAutapase0_iMSNAR0_iMSNAutapase0_FSAR0_FSAutapase1_trial#*.mat');
fidx{2}=dir('Prama_extract_ResultsGD4FS01_Network_isolated0_DA1_dMSNAR0_dMSNAutapase0_iMSNAR0_iMSNAutapase0_FSAR1_FSAutapase1_trial#*.mat');
fidx{3}=dir('Prama_extract_ResultsGD4FS01_Network_isolated0_DA1_dMSNAR0_dMSNAutapase0_iMSNAR0_iMSNAutapase0_FSAR2_FSAutapase1_trial#*.mat');
groupname={'Aut & AR0','Aut & AR1','Aut & AR2'};
paramname={'spikefreq','Burstrate','Intertrain','Burstnumber','Timeinterspike','Burstdura','Burstpercen','thetapower','alfapower','betapower','lowgamapower','highgamapower'};
cellname={'dMSN','iMSN','FS'};
%%
for kk=1:length(fidx)
    
    param{kk}=nan(40,length(paramname),length(cellname));
    
    for k=1:length(fidx{kk})
        filename=fidx{kk}(k).name;
        load(filename)
        disp([num2str(k),'/',num2str(length(fidx)),'--',filename])
        %%
        param{kk}(k,1,:)=spikefreq;
        param{kk}(k,2,:)=Burstrate;
        param{kk}(k,3,:)=Intertrain;
        param{kk}(k,4,:)=Burstnumber;
        param{kk}(k,5,:)=Timeinterspike;
        param{kk}(k,6,:)=Burstdura;
        param{kk}(k,7,:)=Burstpercen;
        param{kk}(k,8,:)=thetapower.*ones(1,3);
        param{kk}(k,9,:)=alfapower.*ones(1,3);
        param{kk}(k,10,:)=betapower.*ones(1,3);
        param{kk}(k,11,:)=lowgamapower.*ones(1,3);
        param{kk}(k,12,:)=highgamapower.*ones(1,3);
        
    end
end

%%

x=1:length(groupname);

for i=1:length(paramname)
    for j=1:length(cellname)
        for kk=1:length(fidx)
            if kk==1
                y=param{kk}(:,i,j);
                ymean=mean(y);
                yn=length(y);
                ysem=std(y)./sqrt(yn);
            elseif kk==2
                Y=param{kk}(:,i,j);
                Ymean=mean(Y);
                Yn=length(Y);
                Ysem=std(Y)./sqrt(Yn);
            elseif kk==3
                YY=param{kk}(:,i,j);
                YYmean=mean(YY);
                YYn=length(YY);
                YYsem=std(YY)./sqrt(YYn);
            end
        end
        figure((j-1)*12+i)
        violinplot([y;Y;YY],[ones(size(y));2*ones(size(Y));3*ones(size(YY))],'ShowData',false);
        ylabel(paramname(i))
        xlabel(cellname(j))
        xlim([0.5 3.5])
        if ymean~=1 & YYmean~=1
            group=[ones(size(y,1),1);2*ones(size(Y,1),1);3*ones(size(YY,1),1)];
            [p((j-1)*12+i),t,stats]=kruskalwallis([y;Y;YY],group,'off');
            c{(j-1)*12+i} = multcompare(stats,'Display','off');
            text(1.5,mean(y),num2str(p((j-1)*12+i)))
        end
    end
end
x=1;
V=[];
X=[];
for i=8:length(paramname)
    for kk=1:length(fidx)
        if kk==1
            y=param{kk}(:,i,1);
            ymean=mean(y);
            yn=length(y);
            ysem=std(y)./sqrt(yn);
        elseif kk==2
            Y=param{kk}(:,i,1);
            Ymean=mean(Y);
            Yn=length(Y);
            Ysem=std(Y)./sqrt(Yn);
        elseif kk==3
            YY=param{kk}(:,i,1);
            YYmean=mean(YY);
            YYn=length(YY);
            YYsem=std(YY)./sqrt(YYn);
        end
    end
    V=[V;y;Y;YY];
    X=[X;x*ones(size(y));(x+1)*ones(size(Y));(x+2)*ones(size(YY))];
    x=x+3;
end
figure(100),clf
violinplot(V,X,'ShowData',false,'Width',0.32);