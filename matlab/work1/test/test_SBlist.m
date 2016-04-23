% purpose:
% to test the correctness of the SB lists.
% polarity_change_dates_different_solar_activity
sector_polarity_change_dates_SBlist
% sector_polarity_change_dates_SBlist_CIR
x=toward_away_date;
y=SBlist(ismember(SBlist(:,2:3),x,'rows'),1);
DailyF107_smooth=DailyF107;
DailyF107_smooth(:,3)=smooth(DailyF107_smooth(:,3),11);
z=DailyF107_smooth(ismember(DailyF107_smooth(:,1:2),x,'rows'),3);

subplot(3,1,1)
plot(x(:,2),'.')
hold on
refline(0,35)
refline(0,125)
refline(0,217)
refline(0,311)
ylabel('doy')
ylim([0,366])
set(gca,'ytick',[82,175,268,357],'yticklabel',['me';'js';'se';'ds'])

subplot(3,1,2)
plot(y,'.')
ylabel('Polarity')
ylim([-1.5,1.5])

subplot(3,1,3)
plot(z,'.')
hold on
ylabel('F10.7')
ylim([50,300])
refline(0,100)
refline(0,160)

xx=unique(x,'rows');
size(x,1)
size(xx,1)

picdir='/home/gdj/study/matlabworkspace/graduation_project/picture/';
picname='test_SBdates.eps';
saveas(gcf,[picdir,picname],'psc2')
