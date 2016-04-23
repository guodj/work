% plot the text with the ratio of the gca
% by Guo, 2015/5/26
function htext = mytext(xr,yr,str,varargin)
tmp = get(gca,'ylim');
ytext = tmp(1) + yr*(tmp(2) - tmp(1));
tmp = get(gca,'xlim');
xtext = tmp(1) + xr*(tmp(2) - tmp(1));

ht = text(xtext,ytext,str,varargin{:});
if nargout,htext=ht; end
end
%-------------------------------------------------------------------------%
% xxtext and yytext is the ratio factor of the gca
%-------------------------------------------------------------------------%
