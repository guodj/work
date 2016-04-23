function sblist_select = is_away_toward( sblist )
%is_away_toward find away to toward SBs
%   input: sblist [year, doy, type]
%   output: logical array
sblist_select = sblist(:,3)==-1;
end

