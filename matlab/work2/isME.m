function medates = isME( dates )
%isME find dates around March equinox
%   input: dates [year, doy, type]
%   output: logical array
%   note that ME:52-142, SE:239-329, JS: 145-235, DS: <53 or >328
%   note that ME:35-125, SE:221-311, JS: 128-218, DS: >311 or <36
medates=dates(:, 2)>35 & dates(:,2)<125;
end
