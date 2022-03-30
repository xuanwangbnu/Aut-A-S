
clear;clc;
fid=dir('dataI_*q200_*b60_*#2.mat');
xlstitle='Result_20060';
freq=200;
num=60;
load('dt.mat');
SRtimewindow=0.004;
stimlatency=0.0015;
prenum=nan(100,1);
posnum=nan(100,1);
%%
for k=1:length(fid)
    filename=fid(k).name;
    xlsname=[filename(1:end-4) '.xlsx'];
    [Data,DataText]=xlsread(xlsname);
    total=Data(2:size(Data,1),2);
    total=total./1000;
    stimtime=1/freq*(num-1)+SRtimewindow;
    pretime=[0,2];
    postime=[2+stimtime,2+stimtime+0.3];
    prenum(k)=length(find(total>pretime(1)&total<pretime(2)));
    posnum(k)=length(find(total>postime(1)&total<postime(2)));
end
 prefreq=nanmean(prenum)/2;
 posfreq=nanmean(posnum)/0.3;
 save([xlstitle,'.mat'])