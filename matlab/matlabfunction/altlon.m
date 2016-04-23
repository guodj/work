function [alt,lon] = altlon( x,y,z )
% calculate altitude and longitude with 
lon=atan2(y,x);
lon=lon/pi*180;
tot=sqrt(x^2+y^2+z^2);
alt=asin(z/tot);
alt=alt/pi*180;


end

