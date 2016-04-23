function lrdates=epoch2lrdates(epoch_date, lrdays)
    % This function returns lrdates around epoch_date
    % input: 
    % epoch_date: such as sblist, CIRlist.
    % lrdays: days before and after epoch date. lrdates=[epochdate-ldays,epochdate+rdays]. Pay attention to the signs of the lrdays.
    % output:
    % lrdates=[left_year, left_doy, right_year, right_doy];
    if isempty(epoch_date)
        lrdates=[];
        return
    end

    if size(lrdays,2)==1
        ldays=lrdays;
        rdays=lrdays;
    elseif size(lrdays,2)==2
        ldays=lrdays(:,1);
        rdays=lrdays(:,2);
    end

    year=epoch_date(:,1);
    doy=epoch_date(:,2);

    left_year=year;
    left_doy=doy-ldays;
    right_year=year;
    right_doy=doy+rdays;
    lyeardoy=numdate(datenum(left_year,1,left_doy));
    ryeardoy=numdate(datenum(right_year,1,right_doy));

    lrdates=[lyeardoy,ryeardoy];
end
