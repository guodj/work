% obtain DailyAp from indexnew
len=size(indexnew,1);
DailyAp=zeros(len,3);
DailyAp(:,1:2)=date2doy(indexnew(:,1:3));
DailyAp(:,3)=nanmean(indexnew(:,4:11),2);