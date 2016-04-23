% obtain percent variation of thermosphere density  during CIR interval
display 'please enter ''yesCIR'' or ''noCIR''!'
sector_polarity_change_dates_SBlist
density_log_actual
% absolute change of mass density
% 250 km
[me_atd250r,me_atd250r_num]=...
    resp2SB(away_toward_date_me,obser250_actual,'perc');
[se_atd250r,se_atd250r_num]=...
    resp2SB(away_toward_date_se,obser250_actual,'perc');
[js_atd250r,js_atd250r_num]=...
    resp2SB(away_toward_date_js,obser250_actual,'perc');
[ds_atd250r,ds_atd250r_num]=...
    resp2SB(away_toward_date_ds,obser250_actual,'perc');
[me_tad250r,me_tad250r_num]=...
    resp2SB(toward_away_date_me,obser250_actual,'perc');
[se_tad250r,se_tad250r_num]=...
    resp2SB(toward_away_date_se,obser250_actual,'perc');
[js_tad250r,js_tad250r_num]=...
    resp2SB(toward_away_date_js,obser250_actual,'perc');
[ds_tad250r,ds_tad250r_num]=...
    resp2SB(toward_away_date_ds,obser250_actual,'perc');
% 400km
[me_atd400r,me_atd400r_num]=...
    resp2SB(away_toward_date_me,obser400_actual,'perc');
[se_atd400r,se_atd400r_num]=...
    resp2SB(away_toward_date_se,obser400_actual,'perc');
[js_atd400r,js_atd400r_num]=...
    resp2SB(away_toward_date_js,obser400_actual,'perc');
[ds_atd400r,ds_atd400r_num]=...
    resp2SB(away_toward_date_ds,obser400_actual,'perc');
[me_tad400r,me_tad400r_num]=...
    resp2SB(toward_away_date_me,obser400_actual,'perc');
[se_tad400r,se_tad400r_num]=...
    resp2SB(toward_away_date_se,obser400_actual,'perc');
[js_tad400r,js_tad400r_num]=...
    resp2SB(toward_away_date_js,obser400_actual,'perc');
[ds_tad400r,ds_tad400r_num]=...
    resp2SB(toward_away_date_ds,obser400_actual,'perc');
% 550km
[me_atd550r,me_atd550r_num]=...
    resp2SB(away_toward_date_me,obser550_actual,'perc');
[se_atd550r,se_atd550r_num]=...
    resp2SB(away_toward_date_se,obser550_actual,'perc');
[js_atd550r,js_atd550r_num]=...
    resp2SB(away_toward_date_js,obser550_actual,'perc');
[ds_atd550r,ds_atd550r_num]=...
    resp2SB(away_toward_date_ds,obser550_actual,'perc');
[me_tad550r,me_tad550r_num]=...
    resp2SB(toward_away_date_me,obser550_actual,'perc');
[se_tad550r,se_tad550r_num]=...
    resp2SB(toward_away_date_se,obser550_actual,'perc');
[js_tad550r,js_tad550r_num]=...
    resp2SB(toward_away_date_js,obser550_actual,'perc');
[ds_tad550r,ds_tad550r_num]=...
    resp2SB(toward_away_date_ds,obser550_actual,'perc');
%%%%%%%%%
figure(3)
yr1=[-30,30];
yr2=[-30,30];
yr3=[-30,30];
h=zeros(1,6);
h(1)=subplot(3,4,1);
myplot_errorbar_text(me_atd250r,me_atd250r_num,...
    me_tad250r,me_tad250r_num,yr1)

h(2)=subplot(3,4,2);
myplot_errorbar_text(se_atd250r,se_atd250r_num,...
    se_tad250r,se_tad250r_num,yr1)

h(3)=subplot(3,4,3);
myplot_errorbar_text(js_atd250r,js_atd250r_num,...
    js_tad250r,js_tad250r_num,yr1)

h(4)=subplot(3,4,4);
myplot_errorbar_text(ds_atd250r,ds_atd250r_num,...
    ds_tad250r,ds_tad250r_num,yr1)

h(5)=subplot(3,4,5);
myplot_errorbar_text(me_atd400r,me_atd400r_num,...
    me_tad400r,me_tad400r_num,yr2)

h(6)=subplot(3,4,6);
myplot_errorbar_text(se_atd400r,se_atd400r_num,...
    se_tad400r,se_tad400r_num,yr2)

h(7)=subplot(3,4,7);
myplot_errorbar_text(js_atd400r,js_atd400r_num,...
    js_tad400r,js_tad400r_num,yr2)

h(8)=subplot(3,4,8);
myplot_errorbar_text(ds_atd400r,ds_atd400r_num,...
    ds_tad400r,ds_tad400r_num,yr2)

h(9)=subplot(3,4,9);
myplot_errorbar_text(me_atd550r,me_atd550r_num,...
    me_tad550r,me_tad550r_num,yr3)

h(10)=subplot(3,4,10);
myplot_errorbar_text(se_atd550r,se_atd550r_num,...
    se_tad550r,se_tad550r_num,yr3)

h(11)=subplot(3,4,11);
myplot_errorbar_text(js_atd550r,js_atd550r_num,...
    js_tad550r,js_tad550r_num,yr3)

h(12)=subplot(3,4,12);
myplot_errorbar_text(ds_atd550r,ds_atd550r_num,...
    ds_tad550r,ds_tad550r_num,yr3)
% title legend xylabel
title(h(1),'March Equinox','fontweight','bold','fontsize',8)
title(h(2),'September Equinox','fontweight','bold','fontsize',8)
title(h(3),'June Solstice','fontweight','bold','fontsize',8)
title(h(4),'December Solstice','fontweight','bold','fontsize',8)
xlabel(h(9),'Epoch Time (Days)','fontweight','bold','fontsize',8);
xlabel(h(10),'Epoch Time (Days)','fontweight','bold','fontsize',8)
xlabel(h(11),'Epoch Time (Days)','fontweight','bold','fontsize',8);
xlabel(h(12),'Epoch Time (Days)','fontweight','bold','fontsize',8)
ylabel(h(1),{'% Change'; '\rho 250 km'},'fontweight','bold','fontsize',8);
ylabel(h(5),{'% Change'; '\rho 400 km'},'fontweight','bold','fontsize',8);
ylabel(h(9),{'% Change'; '\rho 550 km'},'fontweight','bold','fontsize',8);
legend(h(1),'away-toward','toward-away','location','northwest'),
legend(h(1),'boxoff')
% set axes
set(h(1:12),'xtick',-4:2:4)
set(h(1:8),'xticklabel',[]);
set(h([2:4:10,3:4:11,4:4:12]),'yticklabel',[]);
% save picture
set(gcf,'units','centimeter','position',[0,0,21,18])
set(gcf,'units','normalized')
for i=1:12
    a=get (h(i),'position');
    set(h(i),'position',a+[0,-0.04,0.04,0.04])
end
picdir='/home/gdj/study/matlabworkspace/graduation_project/picture/';
picname='relative_thermospheric_density_variation_noCIR.eps';
saveas(gcf,[picdir,picname],'psc2');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
