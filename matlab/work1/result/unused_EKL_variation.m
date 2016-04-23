% merging electric field variation in response to SB crossings.
by=hourlybybzv(:,4);
bz=hourlybybzv(:,5);
v=hourlybybzv(:,6);

bt=sqrt(by.*by+bz.*bz);

bs=bz;bs(bs>0)=0;
bs=abs(bs);

em=[];
em(:,1:3)=hourlybybzv(:,1:3);
em(:,4)=v.*(bt-bz)/(2.*10.^3);
% em unit:mv/m

display 'please enter ''all''!'
sector_polarity_change_dates_SBlist;
% EKL variation
me_athem=hourvalue_change_around_polarity_change(away_toward_date_me,em);
se_athem=hourvalue_change_around_polarity_change(away_toward_date_se,em);
js_athem=hourvalue_change_around_polarity_change(away_toward_date_js,em);
ds_athem=hourvalue_change_around_polarity_change(away_toward_date_ds,em);
me_tahem=hourvalue_change_around_polarity_change(toward_away_date_me,em);
se_tahem=hourvalue_change_around_polarity_change(toward_away_date_se,em);
js_tahem=hourvalue_change_around_polarity_change(toward_away_date_js,em);
ds_tahem=hourvalue_change_around_polarity_change(toward_away_date_ds,em);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% begin drawing picture
figure(1);
hh=[];
hh(1)=subplot(2,1,1);
plot(1:264,me_athem,'b');
hold on;
plot(1:264,me_tahem,'r');

hh(2)=subplot(2,1,2);
plot(1:264,js_athem,'b');
hold on;
plot(1:264,js_tahem,'r');
% set title, label and legend.
legend(hh(1),'away-toward','toward-away','location','northwest'),
legend(hh(1),'boxoff')
title(hh(1),'ME');
ylabel(hh(1),'E_{kl}(mV/m)');
xlabel(hh(1),'Epoch Time (Days)')
title(hh(2),'JS');
ylabel(hh(2),'E_{kl}(mV/m)');
xlabel(hh(2),'Epoch Time (Days)')
% set axes
yr=[0,1.5];
xr=[0,264];
xlim(hh(1),xr),ylim(hh(1),yr)
xlim(hh(2),xr),ylim(hh(2),yr)
set(hh,'xtick',0:24:264,'xticklabel',-5:5)
set(hh,'ytick',0:0.5:1.5)
% save picture
picdir='/home/gdj/study/matlabworkspace/graduation_project/picture/';
picname='hEm.eps';
saveas(gcf,[picdir,picname],'psc2')
