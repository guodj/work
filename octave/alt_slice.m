lat = load('/home/guod/tmp/lat.txt');
lon = load('/home/guod/tmp/lon.txt');
alt = load('/home/guod/tmp/alt.txt');
rho1 = load('/home/guod/tmp/rho1.txt');
rho2 = load('/home/guod/tmp/rho2.txt');

lat = reshape(lat, 72, 36, 50);
lon = reshape(lon, 72, 36, 50);
alt = reshape(alt, 72, 36, 50);
rho1 = reshape(rho1, 72, 36, 50);
rho2 = reshape(rho2, 72, 36, 50);
r = lat+pi/2;
theta = lon;
x = r.*cos(theta);
y = r.*sin(theta);
z = alt;
v = 100*(rho2-rho1)./rho1;
sx = linspace(-pi/2, pi/2, 30);
sy = linspace(-pi/2, pi/2, 30);
sz = z(1,1,30);
slice(x,y,z,rho2,[], [], sz)
