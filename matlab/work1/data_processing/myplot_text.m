function pic = myplot_text( yat,yat_num,yta,yta_num,yr)
%myplot_text(yat,yta,yrl,yru) is just used to simplify the drawing
%process
%   Detailed explanation goes here
xr=-5:5;
yat_prc=prctile(yat,[25,50,75],1);
yta_prc=prctile(yta,[25,50,75],1);
pic=myplot(xr,yat_prc(2,:),'b',xr,yta_prc(2,:),'r');
ylim(yr)
% text
tyup=0.3;
tydown=0.15;
tpy1=tydown*(yr(2)-yr(1))+yr(1);
tpy11=tyup*(yr(2)-yr(1))+yr(1);
text1=num2str(range(yat_prc(2,:)),'%7.1E');
text2=num2str(yat_num);
text3=num2str(range(yta_prc(2,:)),'%7.1E');
text4=num2str(yta_num);
text(-0.5,tpy1,text1,'color','b','fontsize',10);
text(-0.5,tpy11,text3,'color','r','fontsize',10);

end
