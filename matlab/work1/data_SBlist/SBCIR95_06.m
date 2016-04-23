% purpose:
% find all SBs and SBs which do not occur simultaneously with CIRs.
% input: SBlist, streaminterface1
% output:SBlistnoCIR,SBlistyesCIR,SBlist_all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'SBlist67_07.mat'
load 'streaminterface'
dd=1;
% first,change format to datenum
SB_datenum=datenum(SBlist(:,1),1,SBlist(:,2));
stream_datenum=datenum(streaminterface1,'dd-mmm-yyyy');
% second, find CIRs 1 days before and after SBs!!!!!!!!!!!!!!!!!!!!!!
stream_datenum_ex=stream_datenum;
for ii=-dd:dd
    stream_datenum_ex=[stream_datenum_ex;stream_datenum+ii];
end
stream_datenum_ex=unique(stream_datenum_ex);
% third, find SBs without CIRs.
isuncert=SB_datenum<min(stream_datenum_ex) |...
    SB_datenum>max(stream_datenum_ex);
isCIR=ismember(SB_datenum,stream_datenum_ex);
% 3 outputs
SBnoCIR=SBlist(~isuncert & ~isCIR,1:2);
SBall=SBlist(~isuncert,1:2);
SByesCIR=SBlist(~isuncert & isCIR,1:2);
save('F:\mywork\matlabworkspace\ther_dens_sect\data_SBlist\SBtype.mat',...
    'SBnoCIR','SBall','SByesCIR','-append')
