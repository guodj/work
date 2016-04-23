function sblist_select = is_toward_away( sblist )
%ISAT find toward to away SBs
%   input: sblist[year, doy, type]
%   output: logical array
sblist_select = sblist(:,3)==1;
end

