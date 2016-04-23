function pic = plot_median(y1,y2,yr)
%myplot_errorbar_text(yat,yta,yr) is just used to simplify the drawing
%process
xl=5;
xr=-xl:xl;
y11=y1(:,6-xl:6+xl);
y22=y2(:,6-xl:6+xl);
% main plot
y11_prc=prctile(y11,[25,50,75],1);
y22_prc=prctile(y22,[25,50,75],1);
pic=myplot(xr,y11_prc(2,:),'b',xr,y22_prc(2,:),'r');
ylim(gca,yr)
hold on
end


