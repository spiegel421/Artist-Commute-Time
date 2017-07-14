function name = adjust_graphname(name)
% ADJUST_GRAPHNAME Remove suffixes from the graph name

name = name(1:end-5); % truncate .smat
if strendswith(name,'-cc')
    name = name(1:end-3);
elseif strendswith(name,'-wcc')
    name = name(1:end-4);
elseif strendswith(name,'_lcc')
    name = name(1:end-4);
elseif strendswith(name,'-sym')    
    name = name(1:end-4);
end