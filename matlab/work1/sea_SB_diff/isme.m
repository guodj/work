function y = isme( doy )
%isme determines whether doy belongs to ME(+/- 45 days around ME)
me_l=date2doy(1,3,20)-45;
me_r=date2doy(1,3,20)+45;
y=doy>=me_l & doy<=me_r;
end

