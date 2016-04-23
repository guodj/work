% set 999.9=nan, hourly to daily
hour_GSE_Bz(hour_GSE_Bz==999.9)=nan;
daily_GSE_Bz=nan*ones(length(hour_GSE_Bz)/24,3);
hour0=find(hour_GSE_Bz(:,3)==0);
hour23=find(hour_GSE_Bz(:,3)==23);
for irow=1:length(daily_GSE_Bz)
    daily_GSE_Bz(irow,1:2)=hour_GSE_Bz(hour0(irow),1:2);
    daily_GSE_Bz(irow,3)=nanmean(hour_GSE_Bz(hour0(irow):hour23(irow),4));
end
