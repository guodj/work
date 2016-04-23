% set 999.9=nan, hourly to daily
hour_SW_v(hour_SW_v==9999)=nan;
daily_SW_v=nan*ones(length(hour_SW_v)/24,3);
hour0=find(hour_SW_v(:,3)==0);
hour23=find(hour_SW_v(:,3)==23);
for irow=1:length(daily_SW_v)
    daily_SW_v(irow,1:2)=hour_SW_v(hour0(irow),1:2);
    daily_SW_v(irow,3)=nanmean(hour_SW_v(hour0(irow):hour23(irow),4));
end
