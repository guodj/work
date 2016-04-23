function y = isse( doy )
%isse determines whether doy belongs to SE(+/- 45 days around SE)
se_l=date2doy(1,6,20)-45;
se_r=date2doy(1,6,20)+45;
y=doy>=se_l & doy<=se_r;
end

