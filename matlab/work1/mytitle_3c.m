function tit = mytitle_3c( hand,text1,text2,text3,c1,c2,ypos_rat,font_si)
%MYTITLE_3c set title position and default fontsize and fontweight
if nargin==7
    font_si=12;
end
c1_text=num2str(c1);
c2_text=num2str(c2);
tit = title(hand,{[text1,' ( \color[rgb]{',c1_text,'}',text2,...
    '\color[rgb]{0,0,0}, ','\color[rgb]{',c2_text,'}',...
    text3,'\color[rgb]{0,0,0} )']});
yli=get(hand,'ylim');
ypos=yli(1)+ypos_rat*(yli(2)-yli(1));
pp=get(tit,'position');
pp(2)=ypos;
set(tit,'position',pp,'fontweight','bold','fontsize',font_si)


end

