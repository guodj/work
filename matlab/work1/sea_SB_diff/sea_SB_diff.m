% test the seasonal dependence of SB occurrence rate.
load SBlist67_07
% interval:67-77,77-87....
mes=0;ses=0;jss=0;dss=0;
for iyear=1967:10:1997
    isyear=SBlist(:,1)>=iyear & SBlist(:,1)<=iyear+10;
    mesum=sum(isme(SBlist(isyear,2)));
    mes=mes+mesum;
    sesum=sum(isse(SBlist(isyear,2)));
    ses=ses+sesum;
    jssum=sum(isjs(SBlist(isyear,2)));
    jss=jss+jssum;
    dssum=sum(isds(SBlist(isyear,2)));
    dss=dss+dssum;
    xlab=[1,2,3,4];
    y=[mesum,sesum,jssum,dssum];
    bar(xlab,y)
    pause
    
end
% it is found the SE has one more SBs than ME on average.