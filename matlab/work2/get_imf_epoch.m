function imf_epoch=get_imf_epoch(epoch_date,lrdays)
    % get_imf_epoch returns imf values during epoch dates.
    % input:
    % epoch_date=[year,doy]
    % lrdays: days before and after the epoch date
    % output:
    % imf_epoch=[day(from 0),imf];
    % note that the IMF data has 1-hour resolution
    imf_epoch=[];
    if isempty(epoch_date)
        return
    end

    epoch_lrdays=epoch2lrdates(epoch_date,lrdays);
    l=size(epoch_lrdays,1);
    for k=1:l
        imf_temp=get_imf_mdays(epoch_lrdays(k,:));
        imf_epoch=[imf_epoch; imf_temp];
    end
end
