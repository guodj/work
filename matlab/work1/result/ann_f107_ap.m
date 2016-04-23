% study the annual variation of ap
clear
clc
load 'final.mat'
diurnal_avarages_of_f107;
diurnal_averages_of_Ap;
% smooth ap and f107 with time span 365 days
DailyAp_sm=smooth(DailyAp(:,3),365);
DailyF107_sm=smooth(DailyF107(:,3),365);

xdate=DailyAp(:,1:2);
xran=datenum(xdate(:,1),1,xdate(:,2));
% plot
[h,h1,h2]=myplotyy(xran,DailyF107_sm,xran,DailyAp_sm);
% axis set
xticlab=1969:11:2007;
xtic=datenum(xticlab,1,1);
set(h(1),'xtick',xtic,'xticklabel',xticlab)
% title label
ylabel(h(1),'Smoothed F10.7 (10^{-22}W m^{-2} Hz^{-1})')
ylabel(h(2),'Smoothed ap (nT)')
xlabel(h(1),'Year')