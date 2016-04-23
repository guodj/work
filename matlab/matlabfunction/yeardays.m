% calculate the total days in a year.
function y=yeardays(year)
y=year;
y(mod(year,4)~=0)=365;
y(mod(year,100)==0 & mod(year,400)~=0)=365;
y(y~=365)=366;
end