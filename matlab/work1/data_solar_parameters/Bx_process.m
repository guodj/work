% set 999.9=nan, hourly to daily
GSE_Bx(GSE_Bx==999.9)=nan;
daily_GSE_Bx=nan*ones(length(GSE_Bx)/24,3);
hour0=find(GSE_Bx(:,3)==0);
hour23=find(GSE_Bx(:,3)==23);
for irow=1:length(daily_GSE_Bx)
    daily_GSE_Bx(irow,1:2)=GSE_Bx(hour0(irow),1:2);
    daily_GSE_Bx(irow,3)=nanmean(GSE_Bx(hour0(irow):hour23(irow),4));
end
