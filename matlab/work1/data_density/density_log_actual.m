%% obser250,400,550, logarithmetic density to actural density
load density.mat
obser250_actual=obser250(:,1:3);
obser400_actual=obser400(:,1:3);
obser550_actual=obser550(:,1:3);
obser250_actual(:,3)=10.^obser250_actual(:,3)/10^(-12);
obser400_actual(:,3)=10.^obser400_actual(:,3)/10^(-12);
obser550_actual(:,3)=10.^obser550_actual(:,3)/10^(-12);
%% change the time tag of this data set.
obser250_actual(:,3)=circshift(obser250_actual(:,3),[-1,0]);
obser400_actual(:,3)=circshift(obser400_actual(:,3),[-1,0]);
obser550_actual(:,3)=circshift(obser550_actual(:,3),[-1,0]);
%% save density_atual
dir='F:\mywork\matlabworkspace\ther_dens_sect\data_density\';
nam='density.mat';
save([dir,nam],'obser250_actual','obser400_actual',...
    'obser550_actual','-append')