function [ h,h1,h2 ] = myplotyy( x1,y1,x2,y2 )
%MYPLOTYY change some defaults of plotyy
[h,h1,h2]=plotyy(x1,y1,x2,y2);
set(h,'TickLength',[0.01 0.025]);
set(h(2),'xtick',[],'xticklabel',[])
% set(h,'xminortick','on');
% set(h,'yminortick','on');
set([h1,h2],'linewidth',1.5)
set(h(2),'ycolor','r')
set(h2,'color','r')

end

