%-------------------------------------------------------------------------%
% Fix the way that the default dotted line and dash-dot lines look terrible
% when exported as eps
% by Zhong, 2014/11/14
%-------------------------------------------------------------------------%
function fix_dottedline(filename)
fid = fopen(filename,'r');
tempfile = tempname;
outfid = fopen(tempfile,'w');
repeat = 1;
while repeat==1
    thisLine = fgetl(fid);
    iStart = strfind(thisLine,'/DO { [.5');
    if iStart
        thisLine(iStart+7:iStart+8) = '03'; % dotted
    end
    iStart = strfind(thisLine,'/DD { [.5');
    if iStart
        thisLine(iStart+7:iStart+9) = '1.5'; % dotdash 
        thisLine(iStart+10:end+1) = [' ' thisLine(iStart+10:end)];
    end    
    if ~ischar(thisLine)
        repeat = 0;
    else
        fprintf(outfid,'%s\n',thisLine);
    end
end
fclose(fid);
fclose(outfid);
copyfile(tempfile, filename);
delete(tempfile);
end
%-------------------------------------------------------------------------%
% the sample of defination text in *.eps
% line types: solid, dotted, dashed, dotdash
% /SO { [] 0 setdash } bdef
% /DO { [.5 dpi2point mul 4 dpi2point mul] 0 setdash } bdef
% /DA { [6 dpi2point mul] 0 setdash } bdef
% /DD { [.5 dpi2point mul 4 dpi2point mul 6 dpi2point mul 4
%   dpi2point mul] 0 setdash } bdef
%-------------------------------------------------------------------------%
