function [tau, c, d, tr, ts, j, n] = ktau(v, varargin);
% KTAU Compute Kendall's tau and all associated parameters.  
%
% tau = ktau(V) 
% tau = ktau(v1, v2)
%
% Computes the Kendall's tau coefficient between two rankings.  Either v1
% and v2 or V(:,1), V(:,2)
%
% Example:
%  ktau([1 2 3 4],[4 3 2 1],'

%
% David Gleich
% Stanford University
% 8 February 2006
%

%
% Todo
% Update documentation

%
% Updates
% 2008-03-09: Switched to use law-1.3.1 revision and to use truncated
%             precision, updated input parsing to handle truncation case
%

% number of digits to truncate in the computation, -1 uses everything
precision=-1;

firstchararg=find(cellfun(@ischar,varargin),1,'first');
if nargin==1 || (~isempty(firstchararg) && firstchararg == 1),
    v1 = v(:,1); 
    v2 = v(:,2);
else
    v1 = v;
    v2 = varargin{1};
end
if ~isempty(firstchararg), optsu=struct(varargin{firstchararg:end}); else optsu=struct([]); end
if isfield(optsu,'precision') && ~isempty(optsu.precision), precision=optsu.precision; end

try
    tau=compute(v1,v2,precision);
catch
    mypath = which('ktau');
    [pathstr,filename,ext,VERSN] = fileparts(mypath);
    javaaddpath([pathstr filesep 'ktau_lib' filesep 'colt-1.2.0.jar']);
    javaaddpath([pathstr filesep 'ktau_lib' filesep 'fastutil5-5.1.2.jar']);
    javaaddpath([pathstr filesep 'ktau_lib' filesep 'mg4j-2.1.jar']);
    javaaddpath([pathstr filesep 'ktau_lib' filesep 'law-1.3.1.jar']);
    javaaddpath([pathstr filesep 'ktau_lib' filesep 'jsap-2.0.jar']);
    javaaddpath([pathstr filesep 'ktau_lib' filesep 'log4j-1.2.14.jar']);
    javaaddpath([pathstr filesep 'ktau_lib' filesep 'dsiutils-1.0.jar']);    
    
    %res = it.unimi.dsi.law.stat.KendallTau.computeStats(v1, v2, []);
    tau=compute(v1,v2,precision);
end;

% old code to compute tau from raw data
% d = double(res(3));
% tr = double(res(1));
% ts = double(res(2));
% j = double(res(4));
% 
% n = length(v1)*(length(v1)-1)/2;
% 
% c = n - (tr + ts - j) - d;
% 
% tau = (n - (tr + ts - j) - 2*d)/(sqrt(n - tr)*sqrt(n - ts));
end

function rval=compute(v1,v2,precision)
%logger=it.unimi.dsi.Util.getLogger(it.unimi.dsi.law.stat.KendallTau.class);
%logger.
if precision >= 0
    v1 = it.unimi.dsi.law.util.Precision.truncate(v1,precision);
    v2 = it.unimi.dsi.law.util.Precision.truncate(v2,precision);
end
rval = it.unimi.dsi.law.stat.KendallTau.compute(v1, v2);
end