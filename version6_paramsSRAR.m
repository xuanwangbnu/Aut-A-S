
%=====================================================================
% area ---------- cm2
% capacitance --- uF/cm2
% voltage ------- mV
% conductance --- mS/cm2
% current ------- uA/cm2
% time ---------- ms
% concentration --mM

 %% passive properties

Cm=1;

AsFS=1;
AdFS=1;
AdMSN=1;
AiMSN=1;

GsdFS=0.5;% somaDendiCOM

%% channel Conductance 
GLMSN=0.1;
GNaMSN=100;
GKMSN=80;
GMMSN=(DA_flag==0)*1.25+(DA_flag==1)*1.25;

ENaMSN=50;
EKMSN=-100;
ELMSN=-67;

GLFS=0.25;
GNaFS=112.5;
GKFS=225;
GDFS=4;

ENaFS=50;
EKFS=-90;
ELFS=-70;

%% channel dynamic
% MSN
% Na+ current--
alpha_mNaMSN=@(v) nansum([(v~=-54).*0.32.*(v+54)./(1-exp(-(v+54)./4)),(v==-54).*1.28],2);
beta_mNaMSN=@(v) nansum([(v~=-27).*0.28.*(v+27)./(exp((v+27)./5)-1),(v==-27).*1.4],2);
alpha_hNaMSN=@(v) 0.128.*exp(-(v+50)/18);
beta_hNaMSN=@(v) 4./(1+exp(-(v+27)/5));

% Kdr current-- kCurrentMSN
alpha_mKMSN=@(v) nansum([(v~=-52).*0.032.*(v+52)./(1-exp(-(v+52)./5)),(v==-52).*0.16],2);
beta_mKMSN=@(v) 0.5.*exp(-(v+57)/40);

% IM current
Qs=3.2094.*(1e-4);
alpha_mMMSN=@(v) nansum([(v~=-30).*Qs.*(v+30)./(1-exp(-(v+30)/9)),(v==-30).*0.0193],2);
beta_mMMSN=@(v)  nansum([(v~=-30).*-Qs.*(v+30)./(1-exp((v+30)/9)),(v==-30).*2.4538e-05],2);

%% FS
% Na+ current--
inf_mNaFS=@(v) 1./(1+exp(-(v+24)./11.5));
inf_hNaFS=@(v) 1./(1+exp((v+58.3)./6.7));
tau_hNaFS=@(v) 0.5+14./(1+exp((v+60)./12));

% Kdr current-- DB here make neuron stutter
inf_mKFS=@(v) 1./(1+exp(-(v+12.4)./6.8));
tau_mKFS=@(v) (0.087+11.4./(1+exp((v+14.6)./8.6))).*(0.087+11.4./(1+exp(-(v-1.3)./18.7)));

% ID current
inf_mDFS=@(v) 1./(1+exp(-(v+50)/20));
tau_mDFS=2;
inf_hDFS=@(v) 1./(1+exp((v+70)/6));
tau_hDFS=150;

%%
% dendiMultiPoissonExp for FS dendrites
noisyInputPoisson = corrPoisson(Ntypes(3), 50, 2, 1, 1, 2, .5, tmax, dt, 0);
noisyInput_Eesyn = [0];

sigma_noise = 4;% noisyInputMSN

IappFS=(7 + 7.*DA_flag);% dendInput% # 1/5 direct excitation of FSI
IappdMSN=1.19 + 0.1*DA_flag;%injectedCurrentD1 # 4/5 increased excitation to D1 MSNs
IappiMSN=1.19 - 0.1*DA_flag;%injectedCurrentD1 # 5/5 decreased excitation to D2 MSNs

%% synaptic properties
GGAP=0.15 + DA_flag.*(0.3-0.15);%dendDendiGAP，0.33  # 2/5 increased gapjunction conductance
GGABAFS2FS=(0.1 - DA_flag.*(0.1 - 0.005));% somaSomaiSYN，0.58  # 3/5 decreased inhibitory conductance between FSIs
GGABAFSaut=0.1;
GGABAFS2dMSN=0.006;
GGABAFS2iMSN=0.006;
GGABAdMSN2dMSN=0.001;
GGABAiMSN2iMSN=0.001;
EGABA=-80;

if isolated_flag
    GGABAFS2dMSN=0;
    GGABAFS2iMSN=0;
end

%% synaptic properties
% synaptic current parameters
% SynPar{1,1}{1}(1:8)-----  presynaptic parameters
% SynPar{1,1}{2}(1:4)-----  taurise,taudecay
% SynPar{1,1}{3}(1)-------  Synaptic Delay

% ======================= 4~6 mV per synapse,20 synapses per neuron
%                       1   2   3   4      5    6    7  8
% % SynPar{pre,post},% Usr;Tar;Uar;tausr;tauar;taud; NF;x0;
% % ConMtx{post,pre}
SynPar(1,[1,2,3])={{[1,0.01,0.4*dMSNAR_flag,20,150,15,10,5],[0.5,13],round(2/dt)}};% dMSN output 
SynPar(2,[1,2,3])={{[1,0.01,0.4*iMSNAR_flag,20,150,15,10,5],[0.5,13],round(2/dt)}};% iMSN output 
SynPar(3,[1,2,3])    ={{[1,0.01,0.4*0*FSAR_flag,20,150,15,10,5],[0.25,13],round(2/dt)}};% FS output
SynPar(1,4)   =       {{[1,0.01,0.4*2*dMSNAR_flag,20,150,15,10,5],[0.5,13],round(2/dt)}};% iMSN autapse
SynPar(2,4)   =       {{[1,0.01,0.4*2*iMSNAR_flag,20,150,15,10,5],[0.5,13],round(2/dt)}};% dMSN autapse
SynPar(3,4)   =       {{[1,0.01,0.4,20,150*(FSAR_flag+0.01),15,10,5],[0.25,13],round(2/dt)}};% FS autapse

for i=1:length(Ntypes)
    for j=1:length(Ntypes)+1
        %  1   2  3   4   5  6  7   8  9
        % usr,Ca,uar,nsr,nar,N,qsr,qar,r
        SynState(i,j)={{zeros(Ntypes(i),1),...
                        zeros(Ntypes(i),1),...
                        zeros(Ntypes(i),1),...
                        zeros(Ntypes(i),1),...
                        zeros(Ntypes(i),1),...
                        SynPar{i,j}{1}(7).*ones(Ntypes(i),1),...
                        zeros(Ntypes(i),1),...
                        zeros(Ntypes(i),1),...
                        zeros(Ntypes(i),1)}};
        qshift(i,j)={zeros(Ntypes(i),SynPar{i,j}{3})};
        qshift_AR(i,j)={zeros(Ntypes(i),SynPar{i,j}{3})};
        qshift_SR(i,j)={zeros(Ntypes(i),SynPar{i,j}{3})};
    end
end

% residual Ca dynamic.
gamma=1e-3; % amount of per-spike increase in residual Ca,result in ~200 nM/spike increment of Ca++
Ka=800e-3;% affinity of asynchronous release machinery,at least 5 spikes to elicit AR
Kp=800e-3;% affinity of a pump for calcium
beta=20e-3;% maximal rate with which residual calcium is cleared from synaptic terminal by active pumps

Cabaseline=50e-3;
Ip=beta.*Cabaseline.^2./(Cabaseline.^2+Kp.^2);% baseline calcium level
CaO=2e3;

% ---------- RK4 ---------------
stage=4;
coefA=[0,0.5,0.5,1];
coefB=[1/6,2/6,2/6,1/6];

for i=1:length(Ntypes)
    for j=1:length(Ntypes)+1
        KS(i,j)={{zeros(Ntypes(i),stage)}};%sGABA
        S(i,j)={{zeros(Ntypes(i),1)}};%sGABA
        
        KS_AR(i,j)={{zeros(Ntypes(i),stage)}};%sGABA
        S_AR(i,j)={{zeros(Ntypes(i),1)}};%sGABA
        
        KS_SR(i,j)={{zeros(Ntypes(i),stage)}};%sGABA
        S_SR(i,j)={{zeros(Ntypes(i),1)}};%sGABA
    end
end

%% differential equation system
f_mNaMSN=@(mNa,v) alpha_mNaMSN(v).*(1-mNa)-beta_mNaMSN(v).*mNa; 
f_hNaMSN=@(hNa,v) alpha_hNaMSN(v).*(1-hNa)-beta_hNaMSN(v).*hNa; 
f_mKMSN=@(mK,v) alpha_mKMSN(v).*(1-mK)-beta_mKMSN(v).*mK; 
f_mMMSN=@(mM,v) alpha_mMMSN(v).*(1-mM)-beta_mMMSN(v).*mM; 

f_hNaFS=@(hNa,v) 1.*((inf_hNaFS(v)-hNa)./tau_hNaFS(v));
f_mKFS=@(mK,v) 1.*((inf_mKFS(v)-mK)./tau_mKFS(v));
f_mDFS=@(mD,v) (inf_mDFS(v)-mD)./tau_mDFS;
f_hDFS=@(hD,v) (inf_hDFS(v)-hD)./tau_hDFS;

f_vdMSN=@(v,mNa,hNa,mK,mM,sdMSN,sdMSNautapse,sFS,I) ...
    (-AdMSN.*(GLMSN.*(v-ELMSN)...
    +GNaMSN.*mNa.^3.*hNa.*(v-ENaMSN)...
    +GKMSN.*mK.^4.*(v-EKMSN)...
    +GMMSN.*mM.*(v-EKMSN))+I...
    -GGABAdMSN2dMSN*dMSNAutapse_flag.*sdMSNautapse.*(v-EGABA)...
    -GGABAdMSN2dMSN.*ConMtx{1,1}*sdMSN.*(v-EGABA)...
    -GGABAFS2dMSN.*ConMtx{1,3}*sFS.*(v-EGABA))./(AdMSN.*Cm);

f_viMSN=@(v,mNa,hNa,mK,mM,siMSN,siMSNautapse,sFS,I) ...
    (-AiMSN.*(GLMSN.*(v-ELMSN)...
    +GNaMSN.*mNa.^3.*hNa.*(v-ENaMSN)...
    +GKMSN.*mK.^4.*(v-EKMSN)...
    +GMMSN.*mM.*(v-EKMSN))+I...
    -GGABAiMSN2iMSN*iMSNAutapse_flag.*siMSNautapse.*(v-EGABA)...
    -GGABAiMSN2iMSN.*ConMtx{2,2}*siMSN.*(v-EGABA)...
    -GGABAFS2iMSN.*ConMtx{2,3}*sFS.*(v-EGABA))./(AiMSN.*Cm);

f_vsFS=@(vs,vd,hNa,mK,mD,hD,sFS,sFSautapse,mNa,I) ...
    (-AsFS.*(GLFS.*(vs-ELFS)...
    +GNaFS.*mNa.^3.*hNa.*(vs-ENaFS)...
    +GKFS.*mK.^2.*(vs-EKFS)...
    +GDFS.*mD.^3.*hD.*(vs-EKFS))...
    -GsdFS.*(vs-vd)+I...
    -GGABAFSaut*FSAutapse_flag.*sFSautapse.*(vs-EGABA)...
    -GGABAFS2FS.*ConMtx{3,3}*sFS.*(vs-EGABA))./(AsFS.*Cm);

f_vdFS=@(vd,vs,hNa,mK,mD,hD,mNa,I) ...
    (-AdFS.*0.1*(GLFS.*(vd-ELFS)...
    +GNaFS.*mNa.^3.*hNa.*(vd-ENaFS)...
    +GKFS.*mK.^2.*(vd-EKFS)...
    +GDFS.*mD.^3.*hD.*(vd-EKFS))...
    -GsdFS.*(vd-vs)+I...
    -GGAP.*sum(ConMtxGJ{3,3}.*bsxfun(@minus,vd,vd'),2))./(AdFS.*Cm);

f_S=@(S,q,trise,tdecay) -S./tdecay+q.*(1-S)./trise;

%% initial values
vdMSN=-63+63*randn(Ntypes(1),1);%.*ones(Ntypes(1),1);%
mNadMSN=alpha_mNaMSN(vdMSN)./(alpha_mNaMSN(vdMSN)+beta_mNaMSN(vdMSN));
hNadMSN=alpha_hNaMSN(vdMSN)./(alpha_hNaMSN(vdMSN)+beta_hNaMSN(vdMSN));
mKdMSN=alpha_mKMSN(vdMSN)./(alpha_mKMSN(vdMSN)+beta_mKMSN(vdMSN));
mMdMSN=alpha_mMMSN(vdMSN)./(alpha_mMMSN(vdMSN)+beta_mMMSN(vdMSN));
    
viMSN=-63+63*randn(Ntypes(2),1);
mNaiMSN=alpha_mNaMSN(viMSN)./(alpha_mNaMSN(viMSN)+beta_mNaMSN(viMSN));
hNaiMSN=alpha_hNaMSN(viMSN)./(alpha_hNaMSN(viMSN)+beta_hNaMSN(viMSN));
mKiMSN=alpha_mKMSN(viMSN)./(alpha_mKMSN(viMSN)+beta_mKMSN(viMSN));
mMiMSN=alpha_mMMSN(viMSN)./(alpha_mMMSN(viMSN)+beta_mMMSN(viMSN));

vsFS=-90.*ones(Ntypes(3),1);
vdFS=-90.*ones(Ntypes(3),1);

hNa_somaFS=0.05 + 0.85*rand(Ntypes(3),1);
mK_somaFS=0.05 + 0.45*rand(Ntypes(3),1);
mD_somaFS=0.15 + 0.60*rand(Ntypes(3),1);
hD_somaFS=0.15 + 0.45*rand(Ntypes(3),1);

hNa_dendFS=0.05 + 0.85*rand(Ntypes(3),1);
mK_dendFS=0.05 + 0.45*rand(Ntypes(3),1);
mD_dendFS=0.15 + 0.60*rand(Ntypes(3),1);
hD_dendFS=0.15 + 0.45*rand(Ntypes(3),1);

% =====================================================
kmNadMSN=nan(Ntypes(1),stage);
khNadMSN=nan(Ntypes(1),stage);
kmKdMSN=nan(Ntypes(1),stage);
kmMdMSN=nan(Ntypes(1),stage);

kmNaiMSN=nan(Ntypes(2),stage);
khNaiMSN=nan(Ntypes(2),stage);
kmKiMSN=nan(Ntypes(2),stage);
kmMiMSN=nan(Ntypes(2),stage);

khNa_somaFS=nan(Ntypes(3),stage);
kmK_somaFS=nan(Ntypes(3),stage);
kmD_somaFS=nan(Ntypes(3),stage);
khD_somaFS=nan(Ntypes(3),stage);

khNa_dendFS=nan(Ntypes(3),stage);
kmK_dendFS=nan(Ntypes(3),stage);
kmD_dendFS=nan(Ntypes(3),stage);
khD_dendFS=nan(Ntypes(3),stage);

ksdMSN=nan(Ntypes(1),stage);
ksiMSN=nan(Ntypes(2),stage);
ksFS=nan(Ntypes(3),stage);

kvdMSN=nan(Ntypes(1),stage);
kviMSN=nan(Ntypes(2),stage);
kvsFS=nan(Ntypes(3),stage);
kvdFS=nan(Ntypes(3),stage);
% =======================================================
thresh=-25;
spkidxdMSN=[];
spkidxiMSN=[];
spkidxFS=[];

spkidx={spkidxdMSN,spkidxiMSN,spkidxFS};
