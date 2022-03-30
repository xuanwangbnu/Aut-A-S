clear;clc;
fidx{1}=dir('Prama_extract_Results_Network_isolated0_DA1_dMSNAR0_dMSNAutapase0_iMSNAR0_iMSNAutapase0_FSAR0_FSAutapase0_trial#*.mat');
fidx{2}=dir('Prama_extract_Results_Network_isolated0_DA1_dMSNAR0_dMSNAutapase0_iMSNAR0_iMSNAutapase0_FSAR0_FSAutapase1_trial#*.mat');
groupname={'No Aut','Aut & AR'};
paramname={'spikefreq','Burstrate','Intertrain','Burstnumber','Timeinterspike','Burstdura','Burstpercen','thetapower','alfapower','betapower','lowgamapower','highgamapower'};
cellname={'dMSN','iMSN','FS'};
for kk=1:length(fidx)
    param{kk}=nan(10,length(paramname),length(cellname))
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
            end
        end
        figure
        bar(x,[ymean,Ymean],0.24,'FaceColor',[1,1,1],'EdgeColor','k','lineWidth',1),hold on
        errorbar(x,[ymean,Ymean],[ysem,Ysem],'o','color','k','markersize',1,'CapSize',3,'LineWidth',1)
        plot(x-0.1+0.2*rand(size(y,1),1),[y,Y],'ob','MarkerSize',10.5)
        ylabel(paramname(i))
        xlabel(cellname(j))
        
        if  Ymean~=min(Y) & ymean~=min(y);
            if swtest(y,0.05)==0 && swtest(Y,0.05)==0
                [~,p]=ttest2(y,Y);
                title('ttest')
            else
                p=ranksum(y,Y);
                title('ranksum')
            end
            text(1.5,0,num2str(p))
        end
    end
end
x=[1 2];
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
        end
    end
    figure(100)
    bar(x,[ymean,Ymean],0.38,'FaceColor',[1,1,1],'EdgeColor','k','lineWidth',1),hold on
    errorbar(x,[ymean,Ymean],[ysem,Ysem],'o','color','k','markersize',1,'CapSize',3,'LineWidth',1)
    plot(x-0.1+0.2*rand(size(y,1),1),[y,Y],'ob','MarkerSize',12.2)
    ylim([0 0.5])
    x=x+2;
end
