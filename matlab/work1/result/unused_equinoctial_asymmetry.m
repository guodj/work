% study the equinoctial asymmetry of mass density, F10.7 and Ap
% mass density
clear
clc
load 'final.mat'
density_log_actual
ann_density=zeros(3,366);
obs_density=obser400_actual;
for i=1:366
    irow=find(obs_density(:,2)==i);
    ann_density(:,i)=prctile(obs_density(irow,3),[25,50,75]);
end
figure (1)
% ann_density(2,:)=smooth(ann_density(2,:),27);
h1=myplot(1:366,ann_density(2,:),'r');
hold on
h2=myplot(1:366,ann_density(1,:),':');
h3=myplot(1:366,ann_density(3,:),':');
grid on
xlim([1,366])
title('Annual Variation of \rho at 400 km')
xlabel('Doy of Year')
ylabel('\rho, 400 km')
% F10.7
ann_F107=zeros(3,366);
diurnal_avarages_of_f107;
obs_F107=DailyF107;
for i=1:366
    irow=find(obs_F107(:,2)==i);
    ann_F107(:,i)=prctile(obs_F107(irow,3),[25,50,75]);
end
figure (2)
% ann_F107(2,:)=smooth(ann_F107(2,:),27);
h1=myplot(1:366,ann_F107(2,:),'r');
hold on
h2=myplot(1:366,ann_F107(1,:),':');
h3=myplot(1:366,ann_F107(3,:),':');
grid on
xlim([1,366])
title('Annual Variation of F10.7')
xlabel('Doy of Year')
ylabel('F10.7')
% Ap
ann_Ap=zeros(3,366);
diurnal_averages_of_Ap;
obs_Ap=DailyAp;
for i=1:366
    irow=find(obs_Ap(:,2)==i);
    ann_Ap(:,i)=prctile(obs_Ap(irow,3),[25,50,75]);
end
figure (3)
% ann_Ap(2,:)=smooth(ann_Ap(2,:),27);
h1=myplot(1:366,ann_Ap(2,:),'r');
hold on
h2=myplot(1:366,ann_Ap(1,:),':');
h3=myplot(1:366,ann_Ap(3,:),':');
grid on
xlim([1,366])
title('Annual Variation of Ap')
xlabel('Doy of Year')
ylabel('Ap')