% obtain solar parameter variation.
% SBs are clasified as geoei and geoie conditions.
clear
clc
load 'final.mat'
sector_polarity_change_dates_SBlist
diurnal_avarages_of_f107;
diurnal_averages_of_Ap
% variation of F107
% away_toward
SBtype='actual';
% SBtype2='residual';
[f107_ie_eq,f107_ie_eq_num]=resp2SB(SB_ie_eq,DailyF107,SBtype);
[f107_ei_eq,f107_ei_eq_num]=resp2SB(SB_ei_eq,DailyF107,SBtype);
[f107_ie_so,f107_ie_so_num]=resp2SB(SB_ie_so,DailyF107,SBtype);
[f107_ei_so,f107_ei_so_num]=resp2SB(SB_ei_so,DailyF107,SBtype);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Bze,Bzm,Ap responses
[eb_ie_eq,eb_ie_eq_num]=resp2SB(SB_ie_eq,DailyIMFv(:,[1,2,3]),SBtype);
[eb_ei_eq,eb_ei_eq_num]=resp2SB(SB_ei_eq,DailyIMFv(:,[1,2,3]),SBtype);
[eb_ie_so,eb_ie_so_num]=resp2SB(SB_ie_so,DailyIMFv(:,[1,2,3]),SBtype);
[eb_ei_so,eb_ei_so_num]=resp2SB(SB_ei_so,DailyIMFv(:,[1,2,3]),SBtype);
[mb_ie_eq,mb_ie_eq_num]=resp2SB(SB_ie_eq,DailyIMFv(:,[1,2,4]),SBtype);
[mb_ei_eq,mb_ei_eq_num]=resp2SB(SB_ei_eq,DailyIMFv(:,[1,2,4]),SBtype);
[mb_ie_so,mb_ie_so_num]=resp2SB(SB_ie_so,DailyIMFv(:,[1,2,4]),SBtype);
[mb_ei_so,mb_ei_so_num]=resp2SB(SB_ei_so,DailyIMFv(:,[1,2,4]),SBtype);
[ap_ie_eq,ap_ie_eq_num]=resp2SB(SB_ie_eq,DailyAp,SBtype);
[ap_ei_eq,ap_ei_eq_num]=resp2SB(SB_ei_eq,DailyAp,SBtype);
[ap_ie_so,ap_ie_so_num]=resp2SB(SB_ie_so,DailyAp,SBtype);
[ap_ei_so,ap_ei_so_num]=resp2SB(SB_ei_so,DailyAp,SBtype);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% solar wind speed response
[v_ie_eq,v_ie_eq_num]=resp2SB(SB_ie_eq,DailyIMFv(:,[1,2,6]),SBtype);
[v_ei_eq,v_ei_eq_num]=resp2SB(SB_ei_eq,DailyIMFv(:,[1,2,6]),SBtype);
[v_ie_so,v_ie_so_num]=resp2SB(SB_ie_so,DailyIMFv(:,[1,2,6]),SBtype);
[v_ei_so,v_ei_so_num]=resp2SB(SB_ei_so,DailyIMFv(:,[1,2,6]),SBtype);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin to draw pictures
figure(1)
yl1=[80,160];
yl2=[-1.5,1.5];
yl3=[-1.5,1.5];
yl4=[350,550];
yl5=[0,20];
h=[];
h(1)=subplot(5,2,1);
myplot(f107_ie_eq,f107_ei_eq,yl1);
h(2)=subplot(5,2,2);
myplot(f107_ie_so,f107_ei_so,yl1);
h(3)=subplot(5,2,3);
myplot(eb_ie_eq,eb_ei_eq,yl2);
h(4)=subplot(5,2,4);
myplot(eb_ie_so,eb_ei_so,yl2);
h(5)=subplot(5,2,5);
myplot(mb_ie_eq,mb_ei_eq,yl3);
h(6)=subplot(5,2,6);
myplot(mb_ie_so,mb_ei_so,yl3);
h(7)=subplot(5,2,7);
myplot(v_ie_eq,v_ei_eq,yl4);
h(8)=subplot(5,2,8);
myplot(v_ie_so,v_ei_so,yl4);
h(9)=subplot(5,2,9);
myplot(ap_ie_eq,ap_ei_eq,yl5);
h(10)=subplot(5,2,10);
myplot(ap_ie_so,ap_ei_so,yl5);
legend(h(1),'I-E','E-I','location','northwest'),
legend(h(1),'boxoff')
title(h(1),'Equinox','fontweight','bold','fontsize',8)
title(h(2),'Solstice','fontweight','bold','fontsize',8)
set(h(9),'ytick',0:5:15)
set(h(1:10),'xtick',-4:2:4)
set(h(1:8),'xticklabel',[]);
set(h(2:2:10),'yticklabel',[]);
xlabel(h(9),'Epoch Time (Days)','fontweight','bold','fontsize',8);
xlabel(h(10),'Epoch Time (Days)','fontweight','bold','fontsize',8)
ylabel(h(1),'F10.7 (10^{-22}Wm^{-2}Hz^{-1})',...
    'fontweight','bold','fontsize',8);
ylabel(h(3),'GSE Bz (nT)','fontweight','bold','fontsize',8);
ylabel(h(5),'GSM Bz (nT)','fontweight','bold','fontsize',8);
ylabel(h(7),'V_{sw} (km/s)','fontweight','bold','fontsize',8);
ylabel(h(9),'Ap Index (nT)','fontweight','bold','fontsize',8);
% save picture
set(gcf,'units','centimeter','position',[0,0,18,20])
set(gcf,'units','normalized')
for i=1:10
    a=get (h(i),'position');
    set(h(i),'position',a+[0,-0.04,0.03,0.04])
end
% picdir='F:\graduation_project\picture\';
% picname='solar_parameters_variation_geo.eps';
% saveas(gcf,[picdir,picname],'psc2');

