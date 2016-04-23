%-------------------------------------------------------------------------%
% plot the text with the ratio of the gca
% display text as text1(text2,text3).colors of text2 and text3 are
% specified
% by Zhong, 2014/12/3. modified by Guo,2014/12/10
%-------------------------------------------------------------------------%
function htext = plot_text_3c(xxtext,yytext,c1,c2,text1,text2,text3,font_size)
if nargin==7
    font_size=8;
end
tmp = get(gca,'ylim');
ytext = tmp(1) + yytext*(tmp(2) - tmp(1));
tmp = get(gca,'xlim');
xtext = tmp(1) + xxtext*(tmp(2) - tmp(1));
c1_text=num2str(c1);
c2_text=num2str(c2);

htext = text(xtext,ytext,[text1,' ( \color[rgb]{',c1_text,'}',text2,...
    '\color[rgb]{0,0,0}, ','\color[rgb]{',c2_text,'}',...
    text3,'\color[rgb]{0,0,0} )'],...
    'Fontweight','bold','fontsize',font_size);
end
%-------------------------------------------------------------------------%
% xxtext and yytext is the ratio factor of the gca
%-------------------------------------------------------------------------%