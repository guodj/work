function imf_epoch_median=get_imf_epoch_median(imf_epoch)
    % get median values of imf parameters at the same epoch time during some period.
    % input:
    % imf_epoch: imf associated parameters during several epoch period.
    % output:
    % imf_epoch_median=[day(>0),imf associated parameters];=[] when no data exists.
    % note that the IMF data has 1-hour resolution
    imf_epoch=sortrows(imf_epoch,1);
    ldiff=diff(imf_epoch(:,1));
    ldiff=find(ldiff>0.5/24);
    edge=[0;ldiff;size(imf_epoch,1)];

    imf_epoch_median=ones(length(edge)-1,size(imf_epoch,2))*nan;
    for k=1:(length(edge)-1)
        imf_epoch_median(k,:)=nanmedian(imf_epoch((edge(k)+1):edge(k+1),:),1);
    end
end
