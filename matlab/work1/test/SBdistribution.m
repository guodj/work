% purpose:
% test SBs distribution
%%%%%%%%%
%1 annual variation
by=1926;
ey=2014;
annualSB=by:ey;
for iyear=by:ey
    begrow=find(SBlist(:,2)==iyear,1);
    endrow=find(SBlist(:,2)==iyear,1,'last');
    annualSB(2,iyear-by+1)=endrow-begrow+1;
end
plot(annualSB(1,:),annualSB(2,:))
doySBlist=sortrows(SBlist,3);
