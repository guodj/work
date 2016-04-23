%obtain daily percentage change of density 
%dates coverage of gamdmdensitykp2 is larger than obser250
percentdensity=zeros(size(obser250,1),5);
percentdensity(:,1:2)=obser250(:,1:2);
p=find(ismember(gamdmdensitykp2(:,1:2),obser250(:,1:2),'rows'));
percentdensity(:,3)=100*(10.^obser250(:,3)-gamdmdensitykp2(p,3))./gamdmdensitykp2(728:15616,3);
percentdensity(:,4)=100*(10.^obser400(:,3)-gamdmdensitykp2(p,4))./gamdmdensitykp2(728:15616,4);
percentdensity(:,5)=100*(10.^obser550(:,3)-gamdmdensitykp2(p,5))./gamdmdensitykp2(728:15616,5);