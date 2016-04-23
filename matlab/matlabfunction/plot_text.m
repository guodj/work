% plot the text with the ratio of the gca
% by Zhong, 2014/12/3
function htext = plot_text(xxtext,yytext,text_str,font_size)
if nargin==3
    font_size=11;
end
tmp = get(gca,'ylim');
ytext = tmp(1) + yytext*(tmp(2) - tmp(1));
tmp = get(gca,'xlim');
xtext = tmp(1) + xxtext*(tmp(2) - tmp(1));

htext = text(xtext,ytext,text_str,'Fontweight','bold','fontsize',font_size);
end
%-------------------------------------------------------------------------%
% xxtext and yytext is the ratio factor of the gca
%-------------------------------------------------------------------------%