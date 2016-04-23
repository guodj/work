function solar_wind_parameters_epoch=get_solar_wind_parameters_epoch(epoch_date,lrdays)
    % get_solar_wind_parameters_epoch returns solar wind parameters during epoch dates.
    % input:
    % epoch_date=[year,doy, ...]
    % lrdays: days before and after the epoch date
    % output:
    % solar_wind_parameters_epoch=[day(from 0),imf];
    % note that the IMF data has 1-hour resolution
    solar_wind_parameters_epoch=[];
    if isempty(epoch_date)
        return
    end

    epoch_lrdates=epoch2lrdates(epoch_date,lrdays);
    l=size(epoch_lrdates,1);
    for k=1:l
        imf_temp=get_solar_wind_parameters_mdays(epoch_lrdates(k,:));
        solar_wind_parameters_epoch=[solar_wind_parameters_epoch; imf_temp];
    end
end
