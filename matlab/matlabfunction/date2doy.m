function y = date2doy( varargin )
% doy() is used to find the doy of (year,month,day)
if nargin==3
    year=varargin{1};
    month=varargin{2};
    day=varargin{3};
    y=datenum(year,month,day)-datenum(year,1,1)+1;
elseif nargin==1
    dat=varargin{1};
    y=datenum(dat)-datenum(dat(:,1),1,1)+1;  
end

