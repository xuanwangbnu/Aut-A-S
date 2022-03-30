clear;clc;

Ntypes=[100,100,50];
GroupName={'dMSN','iMSN','FS'};

Pconn=xlsread('param_ConMtx.xlsx','Chemical_Pcon');
for i=1:length(Ntypes)% pre
    for j=1:length(Ntypes) % post
        if i==j
            ConMtx(j,i)={genmask(Ntypes(i),Ntypes(j),Pconn(j,i),1,1)};
        else
            ConMtx(j,i)={genmask(Ntypes(i),Ntypes(j),Pconn(j,i),1,0)};
        end
                
    end
end

% ======================================== Gap junction
PconnGJ=xlsread('param_ConMtx.xlsx','Electri_Pcon');
for i=1:length(Ntypes)% pre
    for j=1:length(Ntypes) % post
        if i==j
            ConMtxGJ(j,i)={genmask(Ntypes(i),Ntypes(j),PconnGJ(j,i),0,1)};
        else
            ConMtxGJ(j,i)={genmask(Ntypes(i),Ntypes(j),PconnGJ(j,i),0,0)};
        end
    end
end

% Visualization_Connection(Ntypes,GroupName,ConMtx ,Coor,[],[3],[1,2]);

for i=1:length(Ntypes)
    Coor{i}=rand(Ntypes(i),1);
end

dt=0.01;
tmax=4000;
tspan=0:dt:tmax;


isolated_flag_All=[0];
DA_flag_All=[0];
dMSNAR_flag_All=[0];
dMSNAutapse_flag_All=[0];
iMSNAR_flag_All=[0];
iMSNAutapse_flag_All=[0];
FSAR_flag_All=[1];
FSAutapse_flag_All=[1];
NumTrial=6;

iset=0;

AllRun=who('*_flag_*');
temp=1;
for i=1:length(AllRun)
    temp=temp.*numel(eval(AllRun{i}));
end
AllRun=temp*NumTrial;
%% update states
for TT=1:NumTrial
    for kAA=1:length(isolated_flag_All)
        for kBB=1:length(DA_flag_All)
            for kCC=1:length(dMSNAR_flag_All)
                for kDD=1:length(iMSNAR_flag_All)
                    for kEE=1:length(FSAR_flag_All)
                        for kFF=1:length(dMSNAutapse_flag_All)
                            for kGG=1:length(iMSNAutapse_flag_All)
                                for kHH=1:length(FSAutapse_flag_All)
                                    tic
                                    isolated_flag=isolated_flag_All(kAA);
                                    DA_flag=DA_flag_All(kBB);
                                    dMSNAR_flag=dMSNAR_flag_All(kCC);
                                    iMSNAR_flag=iMSNAR_flag_All(kDD);
                                    FSAR_flag=FSAR_flag_All(kEE);
                                    
                                    dMSNAutapse_flag=dMSNAutapse_flag_All(kFF);
                                    iMSNAutapse_flag=iMSNAutapse_flag_All(kGG);
                                    FSAutapse_flag=FSAutapse_flag_All(kHH);
                                    
                                    version6_paramsSRAR
                                    
                                    iset=iset+1;
                                    
                                    rasterdMSN=[];
                                    rasteriMSN=[];
                                    rasterFS=[];
                                    
                                    for i=1:length(Ntypes)
                                        Recordidx{i}=[1:ceil(Ntypes(i)./10):Ntypes(i)];
                                        vRecord{i}=nan(numel(tspan),numel(Recordidx{i}));
                                    end
                                    
                                    LFPRecord=nan(size(tspan));
                                    LFP_ARSRRecord=nan(numel(tspan),2);% 
                                    vMeanRecord=nan(numel(tspan),length(Ntypes));
                                    
                                    t_display=0;
                                    for i=1:length(tspan)
                                        % recording history values
                                        vdMSN_old=vdMSN;
                                        mNadMSN_old=mNadMSN;
                                        hNadMSN_old=hNadMSN;
                                        mKdMSN_old=mKdMSN;
                                        mMdMSN_old=mMdMSN;
                                        
                                        viMSN_old=viMSN;
                                        mNaiMSN_old=mNaiMSN;
                                        hNaiMSN_old=hNaiMSN;
                                        mKiMSN_old=mKiMSN;
                                        mMiMSN_old=mMiMSN;
                                        
                                        vsFS_old=vsFS;
                                        vdFS_old=vdFS ;
                                        hNaFS_soma_old=hNa_somaFS;
                                        mKFS_soma_old=mK_somaFS;
                                        mDFS_soma_old=mD_somaFS;
                                        hDFS_soma_old=hD_somaFS;
                                        
                                        hNaFS_dend_old=hNa_dendFS;
                                        mKFS_dend_old=mK_dendFS;
                                        mDFS_dend_old=mD_dendFS;
                                        hDFS_dend_old=hD_dendFS;
                                        
                                        SynState_old=SynState;
                                        S_old=S;
                                        S_AR_old=S_AR;% 
                                        S_SR_old=S_SR;% 
                                        
                                        qshift=cellfun(@(x) circshift(x,[0,1]),qshift,'UniformOutput',false);
                                        qshift_AR=cellfun(@(x) circshift(x,[0,1]),qshift_AR,'UniformOutput',false);% 
                                        qshift_SR=cellfun(@(x) circshift(x,[0,1]),qshift_SR,'UniformOutput',false);% 
                                        %-------------update instantaneous gatting variable -------------
                                        mNaFS_soma=inf_mNaFS(vsFS_old);
                                        mNaFS_dend=inf_mNaFS(vdFS_old);
                                        
                                        InoisedMSN=sigma_noise.*randn(Ntypes(1),1).*sqrt(dt);
                                        InoiseiMSN=sigma_noise.*randn(Ntypes(2),1).*sqrt(dt);
                                        InoiseFS=-noisyInputPoisson(:,i).*(vdFS-noisyInput_Eesyn);
                                        
                                        IdMSN=zeros(Ntypes(1),1)+IappdMSN+InoisedMSN;
                                        IiMSN=zeros(Ntypes(2),1)+IappiMSN+InoiseiMSN;
                                        IsomaFS=zeros(Ntypes(3),1);
                                        IdendFS=zeros(Ntypes(3),1)+IappFS+InoiseFS;
                                        
                                        %%
                                        for j=1:length(Ntypes)
                                            if ~isempty(spkidx{j})
                                                
                                                for jj=1:length(Ntypes)+1
                                                    SynState_old{j,jj}{1}(spkidx{j})=SynState_old{j,jj}{1}(spkidx{j})+ SynPar{j,jj}{1}(1).*(1-SynState_old{j,jj}{1}(spkidx{j}));% usr_old
                                                    SynState_old{j,jj}{2}(spkidx{j})=SynState_old{j,jj}{2}(spkidx{j})+ gamma.*log(CaO./SynState_old{j,jj}{2}(spkidx{j}));% Ca_old
                                                    SynState{j,jj}{4}(spkidx{j})=binornd(SynState_old{j,jj}{6}(spkidx{j}),SynState_old{j,jj}{1}(spkidx{j}));% nsr
                                                    SynState_old{j,jj}{6}(spkidx{j})=SynState_old{j,jj}{6}(spkidx{j})-SynState{j,jj}{4}(spkidx{j});% N
                                                    SynState{j,jj}{7}(spkidx{j})=SynPar{j,jj}{1}(8).*SynState{j,jj}{4}(spkidx{j});% qsr
                                                end
                                                
                                            else
                                                
                                                for jj=1:length(Ntypes)+1
                                                    SynState{j,jj}{4}=zeros(Ntypes(j),1);% nsr
                                                    SynState{j,jj}{7}=zeros(Ntypes(j),1);% qsr
                                                end
                                                
                                            end
                                            
                                            
                                            for jj=1:length(Ntypes)+1
                                                SynState{j,jj}{2}=SynState_old{j,jj}{2}+(-beta.*(SynState_old{j,jj}{2}.^2./(SynState_old{j,jj}{2}.^2+Kp^2))+Ip).*dt;% Ca
                                                
                                                temp=max([SynState{j,jj}{2}-Cabaseline,zeros(Ntypes(j),1)],[],2);
                                                SynState_old{j,jj}{3}=SynState_old{j,jj}{3}+SynPar{j,jj}{1}(3).*(1-SynState_old{j,jj}{3}).*(temp.^4./(temp.^4+Ka^4)).*dt;% uar
                                                
                                                SynState{j,jj}{5}=binornd(SynState_old{j,jj}{6},SynState_old{j,jj}{3}.*(rand(size(SynState_old{j,jj}{3}))<SynPar{j,jj}{1}(2)));% nar
                                                SynState{j,jj}{8}=SynPar{j,jj}{1}(8).*SynState{j,jj}{5};% qar
                                                SynState_old{j,jj}{6}=SynState_old{j,jj}{6}-SynState{j,jj}{5};% N
                                                
                                                SynState{j,jj}{9}=binornd(SynPar{j,jj}{1}(7)-SynState_old{j,jj}{6},dt./SynPar{j,jj}{1}(6));% r
                                                
                                                SynState{j,jj}{1}=SynState_old{j,jj}{1}-SynState_old{j,jj}{1}./SynPar{j,jj}{1}(4).*dt;
                                                SynState{j,jj}{3}=SynState_old{j,jj}{3}-SynState_old{j,jj}{3}./SynPar{j,jj}{1}(5).*dt;
                                                SynState{j,jj}{6}=SynState_old{j,jj}{6}+SynState{j,jj}{9};
                                                
                                            end
                                            
                                        end
                                        
                                        %%
                                        for j=1:stage
                                            if j==1
                                                kmNadMSN(:,j)=dt*f_mNaMSN(mNadMSN_old,vdMSN_old);
                                                khNadMSN(:,j)=dt*f_hNaMSN(hNadMSN_old,vdMSN_old);
                                                kmKdMSN(:,j)=dt*f_mKMSN(mKdMSN_old,vdMSN_old);
                                                kmMdMSN(:,j)=dt*f_mMMSN(mMdMSN_old,vdMSN_old);
                                                
                                                kmNaiMSN(:,j)=dt*f_mNaMSN(mNaiMSN_old,viMSN_old);
                                                khNaiMSN(:,j)=dt*f_hNaMSN(hNaiMSN_old,viMSN_old);
                                                kmKiMSN(:,j)=dt*f_mKMSN(mKiMSN_old,viMSN_old);
                                                kmMiMSN(:,j)=dt*f_mMMSN(mMiMSN_old,viMSN_old);
                                                
                                                khNa_somaFS(:,j)=dt*f_hNaFS(hNaFS_soma_old,vsFS_old);
                                                kmK_somaFS(:,j)=dt*f_mKFS(mKFS_soma_old,vsFS_old);
                                                kmD_somaFS(:,j)=dt*f_mDFS(mDFS_soma_old,vsFS_old);
                                                khD_somaFS(:,j)=dt*f_hDFS(hDFS_soma_old,vsFS_old);
                                                
                                                khNa_dendFS(:,j)=dt*f_hNaFS(hNaFS_dend_old,vdFS_old);
                                                kmK_dendFS(:,j)=dt*f_mKFS(mKFS_dend_old,vdFS_old);
                                                kmD_dendFS(:,j)=dt*f_mDFS(mDFS_dend_old,vdFS_old);
                                                khD_dendFS(:,j)=dt*f_hDFS(hDFS_dend_old,vdFS_old);
                                                
                                                kvdMSN(:,j)=dt*f_vdMSN(vdMSN_old,mNadMSN_old,hNadMSN_old,mKdMSN_old,mMdMSN_old,S_old{1,1}{1},S_old{1,4}{1},S_old{3,1}{1},IdMSN);
                                                kviMSN(:,j)=dt*f_viMSN(viMSN_old,mNaiMSN_old,hNaiMSN_old,mKiMSN_old,mMiMSN_old,S_old{2,2}{1},S_old{2,4}{1},S_old{3,2}{1},IiMSN);
                                                kvsFS(:,j)=dt*f_vsFS(vsFS_old,vdFS_old,hNaFS_soma_old,mKFS_soma_old,mDFS_soma_old,hDFS_soma_old,S_old{3,3}{1},S_old{3,4}{1},mNaFS_soma,IsomaFS);
                                                kvdFS(:,j)=dt*f_vdFS(vdFS_old,vsFS_old,hNaFS_dend_old,mKFS_dend_old,mDFS_dend_old,hDFS_dend_old,mNaFS_dend,IdendFS);
                                                
                                                for jj=1:length(Ntypes)
                                                    for jjj=1:length(Ntypes)+1
                                                        for qq=1:length(KS{jj,jjj})
                                                            qshift{jj,jjj}(:,1)=SynState{jj,jjj}{7}+SynState{jj,jjj}{8};% SR + AR
                                                            KS{jj,jjj}{qq}(:,j)=dt*f_S(S_old{jj,jjj}{qq},qshift{jj,jjj}(:,end),SynPar{jj,jjj}{2}(1+2*(qq-1)),SynPar{jj,jjj}{2}(2+2*(qq-1)));
                                                            
                                                            qshift_AR{jj,jjj}(:,1)=SynState{jj,jjj}{8};% AR
                                                            KS_AR{jj,jjj}{qq}(:,j)=dt*f_S(S_AR_old{jj,jjj}{qq},qshift_AR{jj,jjj}(:,end),SynPar{jj,jjj}{2}(1+2*(qq-1)),SynPar{jj,jjj}{2}(2+2*(qq-1))); %                                                            
                                                            
                                                            qshift_SR{jj,jjj}(:,1)=SynState{jj,jjj}{7};% SR
                                                            KS_SR{jj,jjj}{qq}(:,j)=dt*f_S(S_SR_old{jj,jjj}{qq},qshift_SR{jj,jjj}(:,end),SynPar{jj,jjj}{2}(1+2*(qq-1)),SynPar{jj,jjj}{2}(2+2*(qq-1))); %                                                                                                                       
                                                        end
                                                    end
                                                    
                                                end
                                                
                                            else
                                                kmNadMSN(:,j)=dt*f_mNaMSN(mNadMSN_old+kmNadMSN(:,j-1).*coefA(j),vdMSN_old+kvdMSN(:,j-1).*coefA(j));
                                                khNadMSN(:,j)=dt*f_hNaMSN(hNadMSN_old+khNadMSN(:,j-1).*coefA(j),vdMSN_old+kvdMSN(:,j-1).*coefA(j));
                                                kmKdMSN(:,j)=dt*f_mKMSN(mKdMSN_old+kmKdMSN(:,j-1).*coefA(j),vdMSN_old+kvdMSN(:,j-1).*coefA(j));
                                                kmMdMSN(:,j)=dt*f_mMMSN(mMdMSN_old+kmMdMSN(:,j-1).*coefA(j),vdMSN_old+kvdMSN(:,j-1).*coefA(j));
                                                
                                                kmNaiMSN(:,j)=dt*f_mNaMSN(mNaiMSN_old+kmNaiMSN(:,j-1).*coefA(j),viMSN_old+kviMSN(:,j-1).*coefA(j));
                                                khNaiMSN(:,j)=dt*f_hNaMSN(hNaiMSN_old+khNaiMSN(:,j-1).*coefA(j),viMSN_old+kviMSN(:,j-1).*coefA(j));
                                                kmKiMSN(:,j)=dt*f_mKMSN(mKiMSN_old+kmKiMSN(:,j-1).*coefA(j),viMSN_old+kviMSN(:,j-1).*coefA(j));
                                                kmMiMSN(:,j)=dt*f_mMMSN(mMiMSN_old+kmMiMSN(:,j-1).*coefA(j),viMSN_old+kviMSN(:,j-1).*coefA(j));
                                                
                                                khNa_somaFS(:,j)=dt*f_hNaFS(hNaFS_soma_old+khNa_somaFS(:,j-1).*coefA(j),vsFS_old+kvsFS(:,j-1).*coefA(j));
                                                kmK_somaFS(:,j)=dt*f_mKFS(mKFS_soma_old+kmK_somaFS(:,j-1).*coefA(j),vsFS_old+kvsFS(:,j-1).*coefA(j));
                                                kmD_somaFS(:,j)=dt*f_mDFS(mDFS_soma_old+kmD_somaFS(:,j-1).*coefA(j),vsFS_old+kvsFS(:,j-1).*coefA(j));
                                                khD_somaFS(:,j)=dt*f_hDFS(hDFS_soma_old+khD_somaFS(:,j-1).*coefA(j),vsFS_old+kvsFS(:,j-1).*coefA(j));
                                                
                                                khNa_dendFS(:,j)=dt*f_hNaFS(hNaFS_dend_old+khNa_dendFS(:,j-1).*coefA(j),vdFS_old+kvdFS(:,j-1).*coefA(j));
                                                kmK_dendFS(:,j)=dt*f_mKFS(mKFS_dend_old+kmK_dendFS(:,j-1).*coefA(j),vdFS_old+kvdFS(:,j-1).*coefA(j));
                                                kmD_dendFS(:,j)=dt*f_mDFS(mDFS_dend_old+kmD_dendFS(:,j-1).*coefA(j),vdFS_old+kvdFS(:,j-1).*coefA(j));
                                                khD_dendFS(:,j)=dt*f_hDFS(hDFS_dend_old+khD_dendFS(:,j-1).*coefA(j),vdFS_old+kvdFS(:,j-1).*coefA(j));
                                                
                                                
                                                kvdMSN(:,j)=dt*f_vdMSN(vdMSN_old+kvdMSN(:,j-1).*coefA(j),mNadMSN_old+kmNadMSN(:,j-1).*coefA(j),hNadMSN_old+khNadMSN(:,j-1).*coefA(j),mKdMSN_old+kmKdMSN(:,j-1).*coefA(j),mMdMSN_old+kmMdMSN(:,j-1).*coefA(j),S_old{1,1}{1}+KS{1,1}{1}(:,j-1).*coefA(j),S_old{1,4}{1}+KS{1,4}{1}(:,j-1).*coefA(j),S_old{3,1}{1}+KS{3,1}{1}(:,j-1).*coefA(j),IdMSN);
                                                kviMSN(:,j)=dt*f_viMSN(viMSN_old+kviMSN(:,j-1).*coefA(j),mNaiMSN_old+kmNaiMSN(:,j-1).*coefA(j),hNaiMSN_old+khNaiMSN(:,j-1).*coefA(j),mKiMSN_old+kmKiMSN(:,j-1).*coefA(j),mMiMSN_old+kmMiMSN(:,j-1).*coefA(j),S_old{2,2}{1}+KS{2,2}{1}(:,j-1).*coefA(j),S_old{2,4}{1}+KS{2,4}{1}(:,j-1).*coefA(j),S_old{3,2}{1}+KS{3,2}{1}(:,j-1).*coefA(j),IiMSN);
                                                kvsFS(:,j)=dt*f_vsFS(vsFS_old+kvsFS(:,j-1).*coefA(j),vdFS_old+kvdFS(:,j-1).*coefA(j),hNaFS_soma_old+khNa_somaFS(:,j-1).*coefA(j),mKFS_soma_old+kmK_somaFS(:,j-1).*coefA(j),mDFS_soma_old+kmD_somaFS(:,j-1).*coefA(j),hDFS_soma_old+khD_somaFS(:,j-1).*coefA(j),S_old{3,3}{1}+KS{3,3}{1}(:,j-1).*coefA(j),S_old{3,4}{1}+KS{3,4}{1}(:,j-1).*coefA(j),mNaFS_soma,IsomaFS);
                                                kvdFS(:,j)=dt*f_vdFS(vdFS_old+kvdFS(:,j-1).*coefA(j),vsFS_old+kvsFS(:,j-1).*coefA(j),hNaFS_dend_old+khNa_dendFS(:,j-1).*coefA(j),mKFS_dend_old+kmK_dendFS(:,j-1).*coefA(j),mDFS_dend_old+kmD_dendFS(:,j-1).*coefA(j),hDFS_dend_old+khD_dendFS(:,j-1).*coefA(j),mNaFS_dend,IdendFS);
                                                
                                                for jj=1:length(Ntypes)
                                                    for jjj=1:length(Ntypes)+1
                                                        for qq=1:length(KS{jj,jjj})
                                                            qshift{jj,jjj}(:,1)=SynState{jj,jjj}{7}+SynState{jj,jjj}{8};
                                                            KS{jj,jjj}{qq}(:,j)=dt*f_S(S_old{jj,jjj}{qq}+KS{jj,jjj}{qq}(:,j-1).*coefA(j),qshift{jj,jjj}(:,end),SynPar{jj,jjj}{2}(1+2*(qq-1)),SynPar{jj,jjj}{2}(2+2*(qq-1)));
                                                            
                                                            qshift_AR{jj,jjj}(:,1)=SynState{jj,jjj}{8};% 
                                                            KS_AR{jj,jjj}{qq}(:,j)=dt*f_S(S_AR_old{jj,jjj}{qq}+KS_AR{jj,jjj}{qq}(:,j-1).*coefA(j),qshift_AR{jj,jjj}(:,end),SynPar{jj,jjj}{2}(1+2*(qq-1)),SynPar{jj,jjj}{2}(2+2*(qq-1)));
                                                            
                                                            qshift_SR{jj,jjj}(:,1)=SynState{jj,jjj}{7};% 
                                                            KS_SR{jj,jjj}{qq}(:,j)=dt*f_S(S_SR_old{jj,jjj}{qq}+KS_SR{jj,jjj}{qq}(:,j-1).*coefA(j),qshift_SR{jj,jjj}(:,end),SynPar{jj,jjj}{2}(1+2*(qq-1)),SynPar{jj,jjj}{2}(2+2*(qq-1)));
                                                            
                                                        end
                                                    end
                                                end
                                                
                                            end
                                        end
                                        
                                        % 5.Set new solution values
                                        mNadMSN=mNadMSN_old+sum(kmNadMSN.*repmat(coefB,1,1),2);
                                        hNadMSN=hNadMSN_old+sum(khNadMSN.*repmat(coefB,1,1),2);
                                        mKdMSN=mKdMSN_old+sum(kmKdMSN.*repmat(coefB,1,1),2);
                                        mMdMSN=mMdMSN_old+sum(kmMdMSN.*repmat(coefB,1,1),2);
                                        vdMSN=vdMSN_old+sum(kvdMSN.*repmat(coefB,1,1),2);
                                        
                                        mNaiMSN=mNaiMSN_old+sum(kmNaiMSN.*repmat(coefB,1,1),2);
                                        hNaiMSN=hNaiMSN_old+sum(khNaiMSN.*repmat(coefB,1,1),2);
                                        mKiMSN=mKiMSN_old+sum(kmKiMSN.*repmat(coefB,1,1),2);
                                        mMiMSN=mMiMSN_old+sum(kmMiMSN.*repmat(coefB,1,1),2);
                                        viMSN=viMSN_old+sum(kviMSN.*repmat(coefB,1,1),2);
                                        
                                        hNa_somaFS=hNaFS_soma_old+sum(khNa_somaFS.*repmat(coefB,1,1),2);
                                        mK_somaFS=mKFS_soma_old+sum(kmK_somaFS.*repmat(coefB,1,1),2);
                                        mD_somaFS=mDFS_soma_old+sum(kmD_somaFS.*repmat(coefB,1,1),2);
                                        hD_somaFS=hDFS_soma_old+sum(khD_somaFS.*repmat(coefB,1,1),2);
                                        
                                        hNa_dendFS=hNaFS_dend_old+sum(khNa_dendFS.*repmat(coefB,1,1),2);
                                        mK_dendFS=mKFS_dend_old+sum(kmK_dendFS.*repmat(coefB,1,1),2);
                                        mD_dendFS=mDFS_dend_old+sum(kmD_dendFS.*repmat(coefB,1,1),2);
                                        hD_dendFS=hDFS_dend_old+sum(khD_dendFS.*repmat(coefB,1,1),2);
                                        
                                        vsFS=vsFS_old+sum(kvsFS.*repmat(coefB,1,1),2);
                                        vdFS=vdFS_old+sum(kvdFS.*repmat(coefB,1,1),2);
                                        
                                        for j=1:length(Ntypes)
                                            for jj=1:length(Ntypes)+1
                                                for qq=1:length(S{j,jj})
                                                    S{j,jj}{qq}=S_old{j,jj}{qq}+sum(KS{j,jj}{qq}.*repmat(coefB,Ntypes(j),1),2);
                                                    S_AR{j,jj}{qq}=S_AR_old{j,jj}{qq}+sum(KS_AR{j,jj}{qq}.*repmat(coefB,Ntypes(j),1),2);% 
                                                    S_SR{j,jj}{qq}=S_SR_old{j,jj}{qq}+sum(KS_SR{j,jj}{qq}.*repmat(coefB,Ntypes(j),1),2);% 
                                                end
                                            end
                                        end
                                        %% record varible (only this part need to be modifed)
                                        spkidxdMSN=find(vdMSN_old<=thresh&vdMSN>=thresh);
                                        if ~isempty(spkidxdMSN),rasterdMSN=[rasterdMSN;ones(size(spkidxdMSN)).*tspan(i),spkidxdMSN];end
                                        
                                        spkidxiMSN=find(viMSN_old<=thresh&viMSN>=thresh);
                                        if ~isempty(spkidxiMSN),rasteriMSN=[rasteriMSN;ones(size(spkidxiMSN)).*tspan(i),spkidxiMSN];end
                                        
                                        spkidxFS=find(vsFS_old<=thresh&vsFS>=thresh);
                                        if ~isempty(spkidxFS),rasterFS=[rasterFS;ones(size(spkidxFS)).*tspan(i),spkidxFS];end
                                        
                                        spkidx={spkidxdMSN,spkidxiMSN,spkidxFS};
                                        
                                        
                                        LFP= sum([GGABAdMSN2dMSN.*ConMtx{1,1}*S{1,1}{1}.*(vdMSN-EGABA)+GGABAFS2dMSN.*ConMtx{1,3}*S{3,1}{1}.*(vdMSN-EGABA)+GGABAdMSN2dMSN.*S{1,4}{1}*dMSNAutapse_flag.*(vdMSN-EGABA);
                                                  GGABAiMSN2iMSN.*ConMtx{2,2}*S{2,2}{1}.*(viMSN-EGABA)+GGABAFS2iMSN.*ConMtx{2,3}*S{3,2}{1}.*(viMSN-EGABA)+GGABAiMSN2iMSN.*S{2,4}{1}*iMSNAutapse_flag.*(viMSN-EGABA);
                                                  GGABAFS2FS.*ConMtx{3,3}*S{3,3}{1}.*(vsFS-EGABA)+GGABAFSaut.*S{3,4}{1}*FSAutapse_flag.*(vsFS-EGABA)]);
                                        
                                        LFPRecord(i)=LFP;
                                        vMeanRecord(i,:)=[mean(vdMSN),mean(viMSN),mean(vsFS)];
                                        
                                        vRecord{1}(i,:)=vdMSN(Recordidx{1});
                                        vRecord{2}(i,:)=viMSN(Recordidx{2});
                                        vRecord{3}(i,:)=vsFS(Recordidx{3});
                                        
                                        LFP_ARSRRecord(i,:)=sum([GGABAFSaut.*S_SR{3,4}{1}*FSAutapse_flag.*(vsFS-EGABA),GGABAFSaut.*S_AR{3,4}{1}*FSAutapse_flag.*(vsFS-EGABA)]);
                                        
                                        if tspan(i)>=t_display
                                            string=sprintf('%d/%d -- %f percentage...\n',iset,AllRun,tspan(i)/tmax*100);
                                            disp(string)
                                            t_display=t_display+500;
                                        end
                                    end
                                    
                                    raster={rasterdMSN,rasteriMSN,rasterFS};
                                    
                                    save(['Results_Network_isolated',num2str(isolated_flag),...
                                        '_DA',num2str(DA_flag),...
                                        '_dMSNAR',num2str(dMSNAR_flag),...
                                        '_dMSNAutapase',num2str(dMSNAutapse_flag),...
                                        '_iMSNAR',num2str(iMSNAR_flag),...
                                        '_iMSNAutapase',num2str(iMSNAutapse_flag),...
                                        '_FSAR',num2str(FSAR_flag),...
                                        '_FSAutapase',num2str(FSAutapse_flag),...
                                        '_trial#',num2str(TT),'.mat'],'raster','vRecord','vMeanRecord','LFPRecord','LFP_ARSRRecord','tspan','Ntypes',...
                                        'Recordidx','GroupName');
                                        t=toc/60
                                end
                            end
                        end
                    end
                end
            end
        end
    end  
end
