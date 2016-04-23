%in order to calculate GAMDM density, find f107, f107a, doy and kp that
%satisfy the requirements. f107: 3 days before the input date, kp: 12 hours
%prior to the inputdate.
%kp=2.2(average of all the Kp index value).
%fmt of gamdm_index: year,doy,f107,f107a,kp.
len=size(indexnew,1);
gamdm_index=zeros(len,5);
gamdm_index(:,[1,4])=indexnew(:,[1,15]);
gamdm_index(:,3)=circshift(indexnew(:,14),[3,0]);
kp0=reshape(indexnew(:,16:23)',len*8,1)/10.;
kp0=circshift(kp0,[12,0]);
kp0=reshape(kp0,8,len)';
gamdm_index(:,5)=nanmean(kp0,2);
gamdm_index(1,2)=1;
for p=2:len
    gamdm_index(p,2)=gamdm_index(p-1,2)+1;
    if gamdm_index(p,1)==gamdm_index(p-1,1)+1
        gamdm_index(p,2)=1;
    end
end
gamdm_index(isnan(gamdm_index(:,3)),:)=[];
gamdm_index(isnan(gamdm_index(:,4)),:)=[];
gamdm_index(isnan(gamdm_index(:,5)),:)=[];
save '/home/gdj/study/graduation_project/data/processed_data/gamdm_index.txt' -ascii gamdm_index