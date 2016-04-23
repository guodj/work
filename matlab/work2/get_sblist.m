% to get sblist from data file
% input: 
% output: sblist=[ year, doy, sbtype, left_day, right_day ]
%----------------------------------------------------------------%
% file_path='D:\mywork\data\SBlist\';
file_path='/data/SBlist/';
file_name='SBlist.txt';
if ~exist([file_path file_name],'file'), error('Wrong filename input'), end
file_id=fopen([file_path file_name]);
file_data=textscan(file_id,'%s %f %f %f %f %f');
file_close=fclose(file_id);
[str_type, year,month,day,left_day,right_day]=file_data{:};

num_type=zeros(length(str_type),1);
num_type(ismember(str_type,'-,+'))=1;
num_type(ismember(str_type,'+,-'))=-1;

doy=datenum(year,month,day)-datenum(year,1,1)+1;

sblist=[year doy num_type left_day right_day];
