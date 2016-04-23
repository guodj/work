% find independent equinoctial SBs 
% with Dst>-30 for +/- 3 days from SB occurrences.
% indep_SBlist
load 'final.mat' 
% obtain SBlist at F:\mywork\matlabworkspace\ther_dens_sect\data
load 'dst67_07.mat' 
% obtain dst6707 at F:\mywork\matlabworkspace\ther_dens_sect\data
% the dst is hourly averaged
lrinter=3;
indep_SBlist=ones(size(SBlist))*nan;
ik=1;
% dates synchronization
aa=ismember(SBlist(:,2:3),dst6707(:,1:2),'rows');
bb=find(aa,1)+1;
ee=find(aa,1,'last')-1;
SB67_07=SBlist(bb:ee,:);

dst_lt_n30_a=find(dst6707(:,4)<-30);
dst_lt_n30_b=ones(12000000,1)*nan;
for iii=-24*lrinter:24*lrinter
    xxu=1+(iii+24*lrinter)*size(dst_lt_n30_a,1):...
        (iii+24*lrinter+1)*size(dst_lt_n30_a,1);
    dst_lt_n30_b(xxu)=dst_lt_n30_a+iii;
end
dst_lt_n30_b(isnan(dst_lt_n30_b))=[];
dst_lt_n30=unique(dst_lt_n30_b);
dst_gt_n30_7d=setdiff(1:size(dst6707,1),dst_lt_n30);
dst_gt_n30_7d=dst_gt_n30_7d';

dstSB=find(ismember(dst6707(:,1:2),SB67_07(:,2:3),'rows'));
indep_SBlist_a=intersect(dstSB,dst_gt_n30_7d);
indep_SBlist_b=dst6707(indep_SBlist_a,1:2);    
indep_SBlist=unique(indep_SBlist_b,'rows');

index_indp=find(ismember(SBlist(:,2:3),indep_SBlist,'rows'));
indep_SBlist(:,3)=SBlist(index_indp,1);
save('F:\mywork\matlabworkspace\ther_dens_sect\data\indep_SB_67_07.mat',...
    'indep_SBlist')
