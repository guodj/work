% purpose:
% make the used data have the same time interval
%
% variables:
% reduced data:DailyIMFv,gamdmdensity,gamdmdensitykp2,indexnew,hourlybybzv
% referred data:obser250
%
% other instructions:
% all the data about to reduce have longer time intervals than obser250
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
begindate=obser250(1,1:2);
enddate=obser250(length(obser250),1:2);
% DailyIMFv
pb=find(DailyIMFv(:,1)==begindate(1)&DailyIMFv(:,2)==begindate(2));
pe=find(DailyIMFv(:,1)==enddate(1)&DailyIMFv(:,2)==enddate(2));
DailyIMFv=DailyIMFv(pb:pe,:);
% gamdmdensity
pb=find(gamdmdensity(:,1)==begindate(1)&gamdmdensity(:,2)==begindate(2));
pe=find(gamdmdensity(:,1)==enddate(1)&gamdmdensity(:,2)==enddate(2));
gamdmdensity=gamdmdensity(pb:pe,:);
% gamdmdensitykp2
pb=find(gamdmdensitykp2(:,1)==begindate(1)&gamdmdensitykp2(:,2)==begindate(2));
pe=find(gamdmdensitykp2(:,1)==enddate(1)&gamdmdensitykp2(:,2)==enddate(2));
gamdmdensitykp2=gamdmdensitykp2(pb:pe,:);
% hourlybybzv
pb=find(hourlybybzv(:,1)==begindate(1)&hourlybybzv(:,2)==begindate(2),1,'first');
pe=find(hourlybybzv(:,1)==enddate(1)&hourlybybzv(:,2)==enddate(2),1,'last');
hourlybybzv=hourlybybzv(pb:pe,:);
% indexnew
indexdoy=date2doy(indexnew(:,1:3));
pb=find(indexdoy(:,1)==begindate(1)&indexdoy(:,2)==begindate(2));
pe=find(indexdoy(:,1)==enddate(1)&indexdoy(:,2)==enddate(2));
indexnew=indexnew(pb:pe,:);