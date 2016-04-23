% obtain absolute variations of thermosphere density to SB
% SBs are classified into geoei and geoie
clear
clc
load 'final.mat'
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sector_polarity_change_dates_SBlist
density_log_actual
mytype='residual';
% absolute change of mass density
% 250 km
[d250_ie_eq,d250_ie_eq_num]=resp2SB(SB_ie_eq,obser250_actual,mytype);
[d250_ei_eq,d250_ei_eq_num]=resp2SB(SB_ei_eq,obser250_actual,mytype);
[d250_ie_so,d250_ie_so_num]=resp2SB(SB_ie_so,obser250_actual,mytype);
[d250_ei_so,d250_ei_so_num]=resp2SB(SB_ei_so,obser250_actual,mytype);
% 400km
[d400_ie_eq,d400_ie_eq_num]=resp2SB(SB_ie_eq,obser400_actual,mytype);
[d400_ei_eq,d400_ei_eq_num]=resp2SB(SB_ei_eq,obser400_actual,mytype);
[d400_ie_so,d400_ie_so_num]=resp2SB(SB_ie_so,obser400_actual,mytype);
[d400_ei_so,d400_ei_so_num]=resp2SB(SB_ei_so,obser400_actual,mytype);
% 550km
[d550_ie_eq,d550_ie_eq_num]=resp2SB(SB_ie_eq,obser550_actual,mytype);
[d550_ei_eq,d550_ei_eq_num]=resp2SB(SB_ei_eq,obser550_actual,mytype);
[d550_ie_so,d550_ie_so_num]=resp2SB(SB_ie_so,obser550_actual,mytype);
[d550_ei_so,d550_ei_so_num]=resp2SB(SB_ei_so,obser550_actual,mytype);
%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
yr1=[-1,1]*1e-11;
yr2=[-4,4]*1e-13;
yr3=[-4,4]*1e-14;
h=zeros(1,6);
%%%%%%%%%%%%%%
h(1)=subplot(3,2,1);
myplot_text(d250_ie_eq,d250_ie_eq_num,d250_ei_eq,d250_ei_eq_num,yr1)

h(2)=subplot(3,2,2);
myplot_text(d250_ie_so,d250_ie_so_num,d250_ei_so,d250_ei_so_num,yr1)

h(3)=subplot(3,2,3);
myplot_text(d400_ie_eq,d400_ie_eq_num,d400_ei_eq,d400_ei_eq_num,yr2)

h(4)=subplot(3,2,4);
myplot_text(d400_ie_so,d400_ie_so_num,d400_ei_so,d400_ei_so_num,yr2)

h(5)=subplot(3,2,5);
myplot_text(d550_ie_eq,d550_ie_eq_num,d550_ei_eq,d550_ei_eq_num,yr3)

h(6)=subplot(3,2,6);
myplot_text(d550_ie_so,d550_ie_so_num,d550_ei_so,d550_ei_so_num,yr3)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title and xylabel
title(h(1),'Equinox','fontweight','bold','fontsize',8)
title(h(2),'Solstice','fontweight','bold','fontsize',8)
xlabel(h(5),'Epoch Time (Days)','fontweight','bold','fontsize',8);
xlabel(h(6),'Epoch Time (Days)','fontweight','bold','fontsize',8)
ylabel(h(1),'\rho, 250 km (kg/m^3)','fontweight','bold','fontsize',8);
ylabel(h(3),'\rho, 400 km (kg/m^3)','fontweight','bold','fontsize',8);
ylabel(h(5),'\rho, 550 km (kg/m^3)','fontweight','bold','fontsize',8);
legend(h(1),'I-E','E-I','location','northwest'),
legend(h(1),'boxoff')
% xytick
set(h(1:6),'xtick',-4:2:4)
set(h(1:4),'xticklabel',[]);
set(h(2:2:6),'yticklabel',[]);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% positions of the picture and each subplot.
set(gcf,'units','centimeter','position',[0,0,15,18])
set(gcf,'units','normalized')
for i=1:6
    a=get (h(i),'position');
    set(h(i),'position',a+[0,-0.04,0.04,0.04])
end
picdir='F:\graduation_project\picture\';
picname='rho_absolute_change_geo.eps';
saveas(gcf,[picdir,picname],'psc2');
