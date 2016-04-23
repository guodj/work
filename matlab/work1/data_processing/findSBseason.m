function y = findSBseason( x,t )
%SBSeason obtain SB list for season 't'
%   seasons:+/- 40 days around equinoxes and solstices
%   x is the SB list
%   t='me','js','se' or 'ds'

%1,determine the season.
hinterval=40;
switch t
    case 'me'
        pos=x(:,2)>=80-hinterval & x(:,2)<=80+hinterval;
    case 'se'
        pos=x(:,2)>=266-hinterval & x(:,2)<=266+hinterval;
    case 'js'
        pos=x(:,2)>=173-hinterval & x(:,2)<=173+hinterval;
    case 'ds'
        pos=x(:,2)<=356+hinterval-365 | x(:,2)>=356-hinterval;
    otherwise
        error('wrong input')
end
y=x(pos,:);
end

