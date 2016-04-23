% verify indep_SBlist
load 'indep_SB_67_07'
% F:\mywork\matlabworkspace\ther_dens_sect\data\indep_SB_67_07
load 'dst67_07'
% F:\mywork\matlabworkspace\ther_dens_sect\data\indep_SB_67_07.mat
for ii=1:size(indep_SBlist,1)
    year=indep_SBlist(ii,1);
    doy=indep_SBlist(ii,2);
    
    ir=find((dst6707(:,1)==year & dst6707(:,2)==doy),1,'first');
    index_x=ir-24*3:ir+24*3;
    date_x=dst6707(index_x,1:3);
    datenum_x=datenum(date_x(:,1),1,date_x(:,2),date_x(:,3),0,0);
    myplot(datenum_x,dst6707(index_x,4));
    refline(0,-30)
    % set figure
    x_tic_label_a=date_x(:,1:2);
    x_tic_label_b=unique(x_tic_label_a,'rows');
    x_tic_label=x_tic_label_b(:,2);
    x_tic=datenum(x_tic_label_b(:,1),1,x_tic_label_b(:,2));
    set(gca,'xtick',x_tic,'xticklabel',x_tic_label)
    title(gca,['year',num2str(x_tic_label_b(1,1))])
    ylabel(gca,'Dst')
    pause
end