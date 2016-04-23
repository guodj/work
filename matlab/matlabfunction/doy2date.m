function [mm,dd] = doy2date( varargin )
% DOY2DATE converted (year,doy) to (month,day)
% input:
% doy2date(year,doy) or doy2date([year,doy])
% output:
% [mm,dd]:month and day
% written by Guo 15/6/10
if nargin==2
    year=varargin{1};
    doy=varargin{2};
elseif nargin==1
    yeardoy=varargin{1};
    year=yeardoy(:,1);
    doy=yeardoy(:,2);
end
doycyear=[0,31,59,90,120,151,181,212,243,273,304,334,365];
doylyear=[0,31,60,91,121,152,182,213,244,274,305,335,366];
a1=(yeardays(year)==365);
a2=(yeardays(year)==366);
mm=year*nan;
dd=year*nan;
for ii=1:12
    mp1=a1 & doy>doycyear(ii) & doy<=doycyear(ii+1);
    mp2=a2 & doy>doylyear(ii) & doy<=doylyear(ii+1);
    mm(mp1 | mp2)=ii;
    dd(mp1)=doy(mp1)-doycyear(ii);
    dd(mp2)=doy(mp2)-doylyear(ii);
end
end