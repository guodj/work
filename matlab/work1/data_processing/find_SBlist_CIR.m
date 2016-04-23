% purpose:
% find all SBs and SBs which do not occur simultaneously with CIRs.
% input: SBlist, streaminterface1,2,3
% output:SBlistnoCIR,SBlistyesCIR,SBlist_all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load 'final.mat'
load 'streaminterface'
% first,change format to datenum
SB_datenum=datenum(SBlist(:,2),1,SBlist(:,3));
stream_datenum1=datenum(streaminterface1,'dd-mmm-yyyy');
stream_datenum2=datenum(streaminterface2,'dd-mmm-yyyy');
% date interval has to be consistent.combine strint1 and strint2
% densityenddate=datenum(2007,1,279);
% limdate=find(stream_datenum2>max(stream_datenum1) &...
%     stream_datenum2<densityenddate);
% stream_datenum2=stream_datenum2(limdate);
% stream_datenum=[stream_datenum1;stream_datenum2];
stream_datenum=stream_datenum1;
% second, find CIRs 1 days before and after SBs!!!!!!!!!!!!!!!!!!!!!!
dd=1;
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
SBlistnoCIR=SBlist(~isuncert & ~isCIR,:);
SBlist_all=SBlist(~isuncert,:);
SBlistyesCIR=SBlist(~isuncert & isCIR,:);
% some useful information
SBnum=size(SBlist_all,1)
SBallsum=size(stream_datenum,1)
SBnoCIRsum=size(SBlistnoCIR,1)
SByesCIRsum=size(SBlistyesCIR,1)
