% v*Bt variation in response to SB crossings.
by=hourlybybzv(:,4);
bz=hourlybybzv(:,5);
v=hourlybybzv(:,6);
bt=sqrt(by.*by+bz.*bz);
bs=bz;bs(bs>0)=0;
bs=abs(bs);
vbt=[];
vbt(:,1:3)=hourlybybzv(:,1:3);
vbt(:,4)=v.*sqrt(by.*by+bz.*bz)/10.^3;
% vbt unit: mv/m
sector_polarity_change_dates_SBlist;
% vBt variation
me_athvbt=hourvalue_change_around_polarity_change(away_toward_date_me,vbt);
se_athvbt=hourvalue_change_around_polarity_change(away_toward_date_se,vbt);
js_athvbt=hourvalue_change_around_polarity_change(away_toward_date_js,vbt);
ds_athvbt=hourvalue_change_around_polarity_change(away_toward_date_ds,vbt);
me_tahvbt=hourvalue_change_around_polarity_change(toward_away_date_me,vbt);
se_tahvbt=hourvalue_change_around_polarity_change(toward_away_date_se,vbt);
js_tahvbt=hourvalue_change_around_polarity_change(toward_away_date_js,vbt);
ds_tahvbt=hourvalue_change_around_polarity_change(toward_away_date_ds,vbt);
% begin to draw picture
figure(2);
hh=[];
hh(1)=subplot(2,1,1);
plot(1:24*11,me_athvbt,'b');
hold on;
plot(1:24*11,me_tahvbt,'r');

hh(2)=subplot(2,1,2);
plot(1:24*11,js_athvbt,'b');
hold on;
plot(1:24*11,js_tahvbt,'r');
% set title, xylabel, legend. 
legend(hh(1),'away-toward','toward-away','location','northwest'),
legend(hh(1),'boxoff')
title(hh(1),'ME');
ylabel(hh(1),'V*B_T(mv/m)');
xlabel(hh(1),'Epoch Time (Days)')
title(hh(2),'JS');
ylabel(hh(2),'V*B_T(mv/m)');
xlabel(hh(2),'Epoch Time (Days)')
% set axes
yr=[1,2.5];
xr=[0,241];
xlim(hh(1),xr),
ylim(hh(1),yr)
xlim(hh(2),xr),
ylim(hh(2),yr)
set(hh,'xtick',0:24:264,'xticklabel',-5:5)
% save picture
picdir='/home/gdj/study/matlabworkspace/graduation_project/picture/';
picname='hvbt.eps';
saveas(gcf,[picdir,picname],'psc2')