% SB distribution due to different seasons, solar 10.7cm fluxes and
% orientations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
clc
load 'final.mat'
load 'streaminterface.mat'
load SBlist67_07
% 1, determine the SB crossing dates
sector_polarity_change_dates_SBlist;
% 2,obtain the daily average F10.7 (DailyF107) during 19650101~20101231
diurnal_avarages_of_f107;
% 3,smooth f10.7
DailyF107_smooth=DailyF107;
DailyF107_smooth(:,3)=smooth(DailyF107_smooth(:,3),11);
% 4,11-day average F10.7 at SBs
away_toward_date_s=[away_toward_date_me;away_toward_date_se;...
    away_toward_date_js;away_toward_date_ds];
toward_away_date_s=[toward_away_date_me;toward_away_date_se;...
    toward_away_date_js;toward_away_date_ds];
isSBat=ismember(DailyF107_smooth(:,1:2),away_toward_date_s,'rows');
isSBta=ismember(DailyF107_smooth(:,1:2),toward_away_date_s,'rows');
F107_at=DailyF107_smooth(isSBat,:);
F107_ta=DailyF107_smooth(isSBta,:);
SBF107=[F107_at;F107_ta];
% 5,draw pictures.
dateat=datenum(F107_at(:,1),1,F107_at(:,2));
atnum=size(dateat,1);
satnum=num2str(atnum);
dateta=datenum(F107_ta(:,1),1,F107_ta(:,2));
tanum=size(dateta,1);
stanum=num2str(tanum);
% nums of SBs at different solar activity
low_sum=sum(F107_at(:,3)<=100)+sum(F107_ta(:,3)<=100);
med_sum=sum(F107_at(:,3)>100 & F107_at(:,3)<160)+...
    sum(F107_ta(:,3)>100 & F107_ta(:,3)<160);
act_sum=sum(F107_at(:,3)>=160)+sum(F107_ta(:,3)>=160);
strlow_sum=num2str(low_sum);
strmed_sum=num2str(med_sum);
stract_sum=num2str(act_sum);
% nums of SBs at different seasons
me_num=length(findSBseason(SBF107,'me'));
se_num=length(findSBseason(SBF107,'se'));
js_num=length(findSBseason(SBF107,'js'));
ds_num=length(findSBseason(SBF107,'ds'));
% nums of SBs associated with CIR
stream_beg_dat=datenum(streaminterface1(1,:));
len=length(streaminterface1);
stream_end_dat=datenum(streaminterface1(len,:));

% plot
myplot(dateat,F107_at(:,3),'b.')
hold on
myplot(dateta,F107_ta(:,3),'r.')
title('Distribution of SBs')
% 8,set gca
xtl=1967:5:2007;
xt=datenum(xtl,1,1);
set(gca,'xtick',xt,'xticklabel',xtl)
xlabel('Year')
ylabel('11-Day Average F10.7 (10^{-22}Wm^{-2}Hz^{-1})')
l=legend(['A-T SBs(',satnum,')'],['T-A SBs(',stanum,')']);
set(l,'location','northwest')
% add 2 reflines
a=refline(0,100);
b=refline(0,160);
set(a,'color','black','linestyle','-','linewidth',1.5)
set(b,'color','black','linestyle','-','linewidth',1.5)
xlen=xt(length(xt))+200;
ylen1=70;
ylen2=130;
ylen3=190;
text(xlen,ylen1,strlow_sum)
text(xlen,ylen2,strmed_sum)
text(xlen,ylen3,stract_sum)
% add season distribution
xlen1=xt(1);
xlen2=xt(2);
ylen1=260;
ylen2=280;
text(xlen1,ylen1,['JS:',num2str(js_num),';'])
text(xlen2,ylen1,['DS:',num2str(ds_num)])
text(xlen1,ylen2,['ME:',num2str(me_num),';'])
text(xlen2,ylen2,['SE:',num2str(se_num)])
% add CIR occurrence
plot([stream_beg_dat,stream_beg_dat],[50,300],'k--')
plot([stream_end_dat,stream_end_dat],[50,300],'k--')
xlen=xt(7);
ylen1=270;
ylen2=290;
text(xlen,ylen1,'SB_{CIR}:186(~49%)')
text(xlen,ylen2,'SB_{total}:377')
% 9, save picture
dir='F:\graduation_project\picture\';
picname='SB_distribution_ori.eps';
saveas(gcf,[dir,picname],'psc2')

