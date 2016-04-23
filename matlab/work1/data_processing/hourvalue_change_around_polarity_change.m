% calculate hourly value change due to SB crossings.
% x1:polarity change dates,2 columns.
% x2:hourly value. year, doy,hour,value
function y=hourvalue_change_around_polarity_change(x1,x2)
y=[];
x22=x2(1:24:size(x2,1),1:2);
for p1=1:size(x1,1)      
        p2=find(x22(:,1)==x1(p1,1)&x22(:,2)==x1(p1,2));
        if  ~isempty(p2)&&p2+5<=size(x22,1) && p2-5>=1
            y=[y;x2(24*(p2-6)+1:24*(p2+5),4)'];
        end
end
y=nanmedian(y,1);
