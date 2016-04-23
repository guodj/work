function pic = myplot_errorbar_text( yat,yta,yr,texfor)
%myplot_errorbar_text(yat,yta,yr) is just used to simplify the drawing
%process
%   Detailed explanation goes here
xr=-5:5;
% main plot
yat_prc=prctile(yat,[25,50,75],1);
yta_prc=prctile(yta,[25,50,75],1);
pic=myplot(xr,yat_prc(2,:),'b',xr,yta_prc(2,:),'r');
ylim(yr)
hold on
% errorbar
err_tick=30;
ebp=3;
atebup=mean(yat_prc(3,:)-yat_prc(2,:));
ateblow=mean(yat_prc(2,:)-yat_prc(1,:));
h=errorbar(ebp,yat_prc(2,ebp+6),ateblow,atebup,'b','linewidth',1.5);
errorbar_tick(h,err_tick) 

ebp=4;
taebup=mean(yta_prc(3,:)-yta_prc(2,:));
taeblow=mean(yta_prc(2,:)-yta_prc(1,:));
h=errorbar(ebp,yta_prc(2,ebp+6),taeblow,taebup,'r','linewidth',1.5);
errorbar_tick(h,err_tick) 
% text
%text position
tyup=0.3;
tydown=0.15;
tpx1=0;
tpy1=tydown*(yr(2)-yr(1))+yr(1);
tpx2=0;
tpy2=tyup*(yr(2)-yr(1))+yr(1);
%text content
text1=num2str(range(yat_prc(2,:)),texfor);
text3=num2str(range(yta_prc(2,:)),texfor);
text2=num2str(size(yat,1));
text4=num2str(size(yta,1));
text(tpx1,tpy1,text1,'color','b','fontsize',10);
text(tpx2,tpy2,text3,'color','r','fontsize',10);

end
