% purpose: 
% calculate daily average F10.7 at the Earth-Sun distance,save them in
% DailyF107
%%%%%%%%%%%%%%%%%%%%%%%%%%%
DailyF107=zeros(size(indexnew,1),3);
DailyF107(:,3)=indexnew(:,14);
DailyF107(:,1:2)=date2doy(indexnew(:,1:3));