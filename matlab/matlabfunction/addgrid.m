% This function add grid with specified linestyle and color
%---------------------------------------------------------------
function  addgrid( h,color,linestyle )
xt=get(h,'xtick');
yt=get(h,'ytick');
xl=get(h,'xlim');
yl=get(h,'ylim');
xtt=repmat(xt,2,1);
yll=repmat(yl',1,size(xt,2));
ytt=repmat(yt,2,1);
xll=repmat(xl',1,size(yt,2));
hold on
plot(xtt,yll,'color',color,'linestyle',linestyle)
plot(xll,ytt,'color',color,'linestyle',linestyle)

end

