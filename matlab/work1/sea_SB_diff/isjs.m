function y = isjs( doy )
%isjs determines whether doy belongs to JS(+/- 45 days around JS)
js_l=date2doy(1,9,20)-45;
js_r=date2doy(1,9,20)+45;
y=doy>=js_l & doy<=js_r;
end

