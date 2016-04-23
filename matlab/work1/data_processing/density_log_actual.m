% obser250,400,550, logarithmetic density to actural density
obser250_actual=obser250(:,1:3);
obser400_actual=obser400(:,1:3);
obser550_actual=obser550(:,1:3);
obser250_actual(:,3)=10.^obser250_actual(:,3);
obser400_actual(:,3)=10.^obser400_actual(:,3);
obser550_actual(:,3)=10.^obser550_actual(:,3);
% in order to change the time tag of this data set.
obser250_actual(:,3)=circshift(obser250_actual(:,3),[-1,0]);
obser400_actual(:,3)=circshift(obser400_actual(:,3),[-1,0]);
obser550_actual(:,3)=circshift(obser550_actual(:,3),[-1,0]);