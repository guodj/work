% purpose:
% obtain the percent change of thermosphere density during different solar
% activities in ME
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
polarity_change_dates_different_solar_activity
%obtain actual mass density rather than logarithmic density.
density_log_actual
diurnal_avarages_of_f107
diurnal_averages_of_Ap
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%solar active mass density
[se_atd250rh,se_atd250rh_num]=...
    resp2SB(away_toward_date_seh,obser250_actual,'perc');
[se_tad250rh,se_tad250rh_num]=...
    resp2SB(toward_away_date_seh,obser250_actual,'perc');
[se_atd400rh,se_atd400rh_num]=...
    resp2SB(away_toward_date_seh,obser400_actual,'perc');
[se_tad400rh,se_tad400rh_num]=...
    resp2SB(toward_away_date_seh,obser400_actual,'perc');
[se_atd550rh,se_atd550rh_num]=...
    resp2SB(away_toward_date_seh,obser550_actual,'perc');
[se_tad550rh,se_tad550rh_num]=...
    resp2SB(toward_away_date_seh,obser550_actual,'perc');
%solar medium mass density
[se_atd250rm,se_atd250rm_num]=...
    resp2SB(away_toward_date_sem,obser250_actual,'perc');
[se_tad250rm,se_tad250rm_num]=...
    resp2SB(toward_away_date_sem,obser250_actual,'perc');
[se_atd400rm,se_atd400rm_num]=...
    resp2SB(away_toward_date_sem,obser400_actual,'perc');
[se_tad400rm,se_tad400rm_num]=...
    resp2SB(toward_away_date_sem,obser400_actual,'perc');
[se_atd550rm,se_atd550rm_num]=...
    resp2SB(away_toward_date_sem,obser550_actual,'perc');
[se_tad550rm,se_tad550rm_num]=...
    resp2SB(toward_away_date_sem,obser550_actual,'perc');
%solar quiet mass density
[se_atd250rl,se_atd250rl_num]=...
    resp2SB(away_toward_date_sel,obser250_actual,'perc');
[se_tad250rl,se_tad250rl_num]=...
    resp2SB(toward_away_date_sel,obser250_actual,'perc');
[se_atd400rl,se_atd400rl_num]=...
    resp2SB(away_toward_date_sel,obser400_actual,'perc');
[se_tad400rl,se_tad400rl_num]=...
    resp2SB(toward_away_date_sel,obser400_actual,'perc');
[se_atd550rl,se_atd550rl_num]=...
    resp2SB(away_toward_date_sel,obser550_actual,'perc');
[se_tad550rl,se_tad550rl_num]=...
    resp2SB(toward_away_date_sel,obser550_actual,'perc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1=figure(4);
yl1=[-70,70];
yl2=[-70,70];
yl3=[-70,70];
h(1)=subplot(3,3,1);
myplot_errorbar_text(se_atd250rl,se_atd250rl_num,...
    se_tad250rl,se_tad250rl_num,yl1)

h(2)=subplot(3,3,2);
myplot_errorbar_text(se_atd250rm,se_atd250rm_num,...
    se_tad250rm,se_tad250rm_num,yl1)

h(3)=subplot(3,3,3);
myplot_errorbar_text(se_atd250rh,se_atd250rh_num,...
    se_tad250rh,se_tad250rh_num,yl1)

h(4)=subplot(3,3,4);
myplot_errorbar_text(se_atd400rl,se_atd400rl_num,...
    se_tad400rl,se_tad400rl_num,yl2)

h(5)=subplot(3,3,5);
myplot_errorbar_text(se_atd400rm,se_atd400rm_num,...
    se_tad400rm,se_tad400rm_num,yl2)

h(6)=subplot(3,3,6);
myplot_errorbar_text(se_atd400rh,se_atd400rh_num,...
    se_tad400rh,se_tad400rh_num,yl2)

h(7)=subplot(3,3,7);
myplot_errorbar_text(se_atd550rl,se_atd550rl_num,...
    se_tad550rl,se_tad550rl_num,yl3)

h(8)=subplot(3,3,8);
myplot_errorbar_text(se_atd550rm,se_atd550rm_num,...
    se_tad550rm,se_tad550rm_num,yl3)

h(9)=subplot(3,3,9);
myplot_errorbar_text(se_atd550rh,se_atd550rh_num,...
    se_tad550rh,se_tad550rh_num,yl3)

% legend title xylabel
legend(h(1),'away-toward','toward-away','location','northwest'),
legend(h(1),'boxoff')
title(h(1),'Solar Minimum','fontweight','bold','fontsize',8)
title(h(2),{'Setpember Equinox';'Solar Medium'},...
    'fontweight','bold','fontsize',8)
title(h(3),'Solar Maximum','fontweight','bold','fontsize',8)
xlabel(h(7),'Epoch Time(Days)','fontweight','bold','fontsize',8);
xlabel(h(8),'Epoch Time(Days)','fontweight','bold','fontsize',8)
xlabel(h(9),'Epoch Time(Days)','fontweight','bold','fontsize',8);
ylabel(h(1),{'%Change'; '\rho 250km'},'fontweight','bold','fontsize',8);
ylabel(h(4),{'%Change'; '\rho 400km'},'fontweight','bold','fontsize',8);
ylabel(h(7),{'%Change'; '\rho 550km'},'fontweight','bold','fontsize',8);
% axes
set(h(1:9),'xtick',-4:2:4)
set(h(1:6),'xticklabel',[]);
set(h([2:3:8,3:3:9]),'yticklabel',[]);
% save picture
set(gcf,'units','centimeter','position',[0,0,18,15])
set(gcf,'units','normalized')
for i=1:9
    a=get (h(i),'position');
    set(h(i),'position',a+[0,-0.05,0.05,0.06])
end
picdir='/home/gdj/study/matlabworkspace/graduation_project/picture/';
picname='relative_thermospheric_density_variation_sesolar.eps';
saveas(gcf,[picdir,picname],'psc2');
