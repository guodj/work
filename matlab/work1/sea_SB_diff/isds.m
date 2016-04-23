function y = isds( doy )
%isds determines whether doy belongs to DS(+/- 45 days around DS)
ds_l=date2doy(1,12,20)-45;
ds_r=date2doy(1,12,20)+45-365;
y=doy>=ds_l | doy<=ds_r;
end

