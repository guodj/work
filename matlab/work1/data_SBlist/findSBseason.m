function y = findSBseason( x,t )
%SBSeason obtain SB list for season 't'
%   seasons:+/- 40 days around equinoxes and solstices
%   x is the SB list
%   t='me','js','se' or 'ds'

%1,determine the season.
hinterval=45;
switch t
    case 'me'
        pos=x(:,2)>=79-hinterval & x(:,2)<=79+hinterval+1;
    case 'se'
        pos=x(:,2)>=263-hinterval-1 & x(:,2)<=263+hinterval;
    case 'js'
        pos=x(:,2)>=171-hinterval & x(:,2)<=171+hinterval;
    case 'ds'
        pos=x(:,2)<=354+hinterval-365-1 | x(:,2)>=354-hinterval;
    otherwise
        error('wrong input')
end
y=x(pos,:);
end

