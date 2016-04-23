% obtain percent variation of thermosphere density
% SBs are classified into geoei and geoie
sector_polarity_change_dates_SBlist
density_log_actual
mytype='perc';
% absolute change of mass density
% 250 km
[d250r_ie_eq,d250r_ie_eq_num]=resp2SB(SB_ie_eq,obser250_actual,mytype);
[d250r_ei_eq,d250r_ei_eq_num]=resp2SB(SB_ei_eq,obser250_actual,mytype);
[d250r_ie_so,d250r_ie_so_num]=resp2SB(SB_ie_so,obser250_actual,mytype);
[d250r_ei_so,d250r_ei_so_num]=resp2SB(SB_ei_so,obser250_actual,mytype);
% 400km
[d400r_ie_eq,d400r_ie_eq_num]=resp2SB(SB_ie_eq,obser400_actual,mytype);
[d400r_ei_eq,d400r_ei_eq_num]=resp2SB(SB_ei_eq,obser400_actual,mytype);
[d400r_ie_so,d400r_ie_so_num]=resp2SB(SB_ie_so,obser400_actual,mytype);
[d400r_ei_so,d400r_ei_so_num]=resp2SB(SB_ei_so,obser400_actual,mytype);
% 550km
[d550r_ie_eq,d550r_ie_eq_num]=resp2SB(SB_ie_eq,obser550_actual,mytype);
[d550r_ei_eq,d550r_ei_eq_num]=resp2SB(SB_ei_eq,obser550_actual,mytype);
[d550r_ie_so,d550r_ie_so_num]=resp2SB(SB_ie_so,obser550_actual,mytype);
[d550r_ei_so,d550r_ei_so_num]=resp2SB(SB_ei_so,obser550_actual,mytype);
%%%%%%%%%%%%%%%%%%%%%%%
figure(3)
yr1=[-30,30];
yr2=[-30,30];
yr3=[-30,30];
h=zeros(1,6);
%%%%%%%%%%%%%%
h(1)=subplot(3,2,1);
myplot_errorbar_text(d250r_ie_eq,d250r_ie_eq_num,d250r_ei_eq,d250r_ei_eq_num,yr1)

h(2)=subplot(3,2,2);
myplot_errorbar_text(d250r_ie_so,d250r_ie_so_num,d250r_ei_so,d250r_ei_so_num,yr1)

h(3)=subplot(3,2,3);
myplot_errorbar_text(d400r_ie_eq,d400r_ie_eq_num,d400r_ei_eq,d400r_ei_eq_num,yr2)

h(4)=subplot(3,2,4);
myplot_errorbar_text(d400r_ie_so,d400r_ie_so_num,d400r_ei_so,d400r_ei_so_num,yr2)

h(5)=subplot(3,2,5);
myplot_errorbar_text(d550r_ie_eq,d550r_ie_eq_num,d550r_ei_eq,d550r_ei_eq_num,yr3)

h(6)=subplot(3,2,6);
myplot_errorbar_text(d550r_ie_so,d550r_ie_so_num,d550r_ei_so,d550r_ei_so_num,yr3)

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% title and xylabel
title(h(1),'Equinox','fontweight','bold','fontsize',8)
title(h(2),'Solstice','fontweight','bold','fontsize',8)
xlabel(h(5),'Epoch Time (Days)','fontweight','bold','fontsize',8);
xlabel(h(6),'Epoch Time (Days)','fontweight','bold','fontsize',8)
ylabel(h(1),{'% Change';'\rho, 250 km'},'fontweight','bold','fontsize',8);
ylabel(h(3),{'% Change';'\rho, 400 km'},'fontweight','bold','fontsize',8);
ylabel(h(5),{'% Change';'\rho, 550 km'},'fontweight','bold','fontsize',8);
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
picname='rho_percent_change_geo.eps';
saveas(gcf,[picdir,picname],'psc2');
