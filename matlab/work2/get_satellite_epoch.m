function [satellite_epoch,varargout]=get_satellite_epoch(epoch_date,lrdays,satellite)
    % get_rho_epoch returns champ, grace or goce density during epoch dates.
    % input:
    % epoch_date=[year,doy]
    % lrdays: days before and after the epoch date
    % satellite='champ', 'grace' or 'goce'
    % output:
    % satellite_epoch: same output as get_champ_grace_mdays or get_goce_mdays.
    satellite_epoch=[];
    if nargout==2
        varargout{1}=[];
    end
    if isempty(epoch_date)
        return
    end

    epoch_lrdays=epoch2lrdates(epoch_date,lrdays);
    l=size(epoch_lrdays,1);
    for k=1:l
        switch satellite
            case {'champ','grace'}
                [uporbit,downorbit]=get_champ_grace_mdays(epoch_lrdays(k,:),satellite); 
            case 'goce'
                [uporbit,downorbit]=get_goce_mdays(epoch_lrdays(k,:));
            otherwise
                error('wrong input of satellite')
        end
        if isempty(uporbit) | isempty(downorbit)
            continue
        end
        satellite_epoch=[satellite_epoch;uporbit;downorbit];

        if nargout==2
            switch satellite
                case {'champ', 'grace'}
                    uprlt=get_rltrho(uporbit(:,[4,13]));
                    downrlt=get_rltrho(downorbit(:,[4,13]));
                case 'goce'
                    uprlt=get_rltrho(uporbit(:,[4,6]));
                    downrlt=get_rltrho(downorbit(:,[4,6]));
                otherwise
                    error('wrong input of satellite')
            end
            varargout{1}=[varargout{1};uprlt;downrlt];
        end
    end
end

function rtlrho=get_rltrho(rho_mdays)
    % rho_mdays=[lat,rho]
    % rtlrho=[mean density at the same latitude, relative_density]
    % note that data in rho_mdays has the same local time
    rtlrho=[];
    if isempty(rho_mdays)
        return
    end

    l=size(rho_mdays,1);
    rtlrho=ones(l,2)*nan;
    s=unique(rho_mdays(:,1));
    for k=1:length(s)
        index = (rho_mdays(:,1)==s(k));
        if sum(index)>l/length(s)/5
            rtlrho(index,1)=mean(rho_mdays(index,2));
            rtlrho(index,2)=((rho_mdays(index,2)-rtlrho(index,1))./rtlrho(index,1))*100;
        end
    end
end
