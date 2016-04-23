% set 999.9=nan, hourly to daily
hour_GSM_Bz(hour_GSM_Bz==999.9)=nan;
daily_GSM_Bz=nan*ones(length(hour_GSM_Bz)/24,3);
hour0=find(hour_GSM_Bz(:,3)==0);
hour23=find(hour_GSM_Bz(:,3)==23);
for irow=1:length(daily_GSM_Bz)
    daily_GSM_Bz(irow,1:2)=hour_GSM_Bz(hour0(irow),1:2);
    daily_GSM_Bz(irow,3)=nanmean(hour_GSM_Bz(hour0(irow):hour23(irow),4));
end
