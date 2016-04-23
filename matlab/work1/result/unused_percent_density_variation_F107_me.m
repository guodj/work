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
[me_atd250rh,me_atd250rh_num]=...
    resp2SB(away_toward_date_meh,obser250_actual,'perc');
[me_tad250rh,me_tad250rh_num]=...
    resp2SB(toward_away_date_meh,obser250_actual,'perc');
[me_atd400rh,me_atd400rh_num]=...
    resp2SB(away_toward_date_meh,obser400_actual,'perc');
[me_tad400rh,me_tad400rh_num]=...
    resp2SB(toward_away_date_meh,obser400_actual,'perc');
[me_atd550rh,me_atd550rh_num]=...
    resp2SB(away_toward_date_meh,obser550_actual,'perc');
[me_tad550rh,me_tad550rh_num]=...
    resp2SB(toward_away_date_meh,obser550_actual,'perc');
%solar medium mass density
[me_atd250rm,me_atd250rm_num]=...
    resp2SB(away_toward_date_mem,obser250_actual,'perc');
[me_tad250rm,me_tad250rm_num]=...
    resp2SB(toward_away_date_mem,obser250_actual,'perc');
[me_atd400rm,me_atd400rm_num]=...
    resp2SB(away_toward_date_mem,obser400_actual,'perc');
[me_tad400rm,me_tad400rm_num]=...
    resp2SB(toward_away_date_mem,obser400_actual,'perc');
[me_atd550rm,me_atd550rm_num]=...
    resp2SB(away_toward_date_mem,obser550_actual,'perc');
[me_tad550rm,me_tad550rm_num]=...
    resp2SB(toward_away_date_mem,obser550_actual,'perc');
%solar quiet mass density
[me_atd250rl,me_atd250rl_num]=...
    resp2SB(away_toward_date_mel,obser250_actual,'perc');
[me_tad250rl,me_tad250rl_num]=...
    resp2SB(toward_away_date_mel,obser250_actual,'perc');
[me_atd400rl,me_atd400rl_num]=...
    resp2SB(away_toward_date_mel,obser400_actual,'perc');
[me_tad400rl,me_tad400rl_num]=...
    resp2SB(toward_away_date_mel,obser400_actual,'perc');
[me_atd550rl,me_atd550rl_num]=...
    resp2SB(away_toward_date_mel,obser550_actual,'perc');
[me_tad550rl,me_tad550rl_num]=...
    resp2SB(toward_away_date_mel,obser550_actual,'perc');
%%%%%%%%%%%%%%%%%%%%%%%%%%%
f1=figure(4);
yl1=[-20,20];
yl2=[-20,20];
yl3=[-20,20];
h(1)=subplot(3,3,1);
myplot_errorbar_text(me_atd250rl,me_atd250rl_num,...
    me_tad250rl,me_tad250rl_num,yl1)

h(2)=subplot(3,3,2);
myplot_errorbar_text(me_atd250rm,me_atd250rm_num,...
    me_tad250rm,me_tad250rm_num,yl1)

h(3)=subplot(3,3,3);
myplot_errorbar_text(me_atd250rh,me_atd250rh_num,...
    me_tad250rh,me_tad250rh_num,yl1)

h(4)=subplot(3,3,4);
myplot_errorbar_text(me_atd400rl,me_atd400rl_num,...
    me_tad400rl,me_tad400rl_num,yl2)

h(5)=subplot(3,3,5);
myplot_errorbar_text(me_atd400rm,me_atd400rm_num,...
    me_tad400rm,me_tad400rm_num,yl2)

h(6)=subplot(3,3,6);
myplot_errorbar_text(me_atd400rh,me_atd400rh_num,...
    me_tad400rh,me_tad400rh_num,yl2)

h(7)=subplot(3,3,7);
myplot_errorbar_text(me_atd550rl,me_atd550rl_num,...
    me_tad550rl,me_tad550rl_num,yl3)

h(8)=subplot(3,3,8);
myplot_errorbar_text(me_atd550rm,me_atd550rm_num,...
    me_tad550rm,me_tad550rm_num,yl3)

h(9)=subplot(3,3,9);
myplot_errorbar_text(me_atd550rh,me_atd550rh_num,...
    me_tad550rh,me_tad550rh_num,yl3)

% legend title xylabel
legend(h(1),'away-toward','toward-away','location','northwest'),
legend(h(1),'boxoff')
title(h(1),'Solar Minimum','fontweight','bold','fontsize',8)
title(h(2),{'March Equinox';'Solar Medium'},'fontweight','bold','fontsize',8)
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
picname='relative_thermospheric_density_variation_mesolar.eps';
saveas(gcf,[picdir,picname],'psc2');
