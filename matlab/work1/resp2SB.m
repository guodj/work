% calculate the changes of a physical value in response to SB crossings in
% different seasons.
% x1: SBlist. 2 columns: year doy
% x2: physical value (mass density, f107, ap...), 3 columns: year, doy,
% value.
% x3:'actual', 'residual', 'perc'
%------------------------------------------------------------------------%
function y=resp2SB(x1,x2,x3)
y=[];
% find positions of SB crossing dates in x2
p=find(ismember(x2(:,1:2),x1,'rows'));
p=p';
time_day=5;
p(2,:)=p(1,:)-time_day;
p(3,:)=p(1,:)+time_day;
% delete columns exceeding matrix dimensions.
dc=p(2,:)<1 | p(3,:)>size(x2,1);
p(:,dc)=[];
%-----------------------------------------------------------------------%
for i=1:size(p,2)
    yadd=x2(p(2,i):p(3,i),3)';
    y=[y;yadd];
end
switch x3
    case 'actual'
    case 'residual'
        y=y-repmat(nanmean(y,2),1,size(y,2));
    case 'perc'
        rowmean=repmat(nanmean(y,2),1,size(y,2));
        y=100*(y-rowmean)./rowmean;
    otherwise
        error('wrong inpur!')
end
