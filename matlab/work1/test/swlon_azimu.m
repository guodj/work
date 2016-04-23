% test whether solar wind longitude is the solar wind azimuthal angle
dn_swlon=zeros(length(solarwindlongitude),1);
dn_swlon(1)=datenum('1964 1 0','yyyy dd hh');
dn_001=datenum('0 0 1','yyyy dd hh');
for irow=2:length(solarwindlongitude)
    dn_swlon
