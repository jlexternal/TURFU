function [seqcat] = get_seqcat(taskid,epimap,seqdir,resp)
%  GET_SEQCAT  Get sequence category
%
%  Usage: [seqcat] = GET_SEQCAT(taskid,epimap,seqdir,[resp])
%
%  where taskid - task identifier => 1:observer or 2:actor
%        epimap - episode mapping => 1:orange=left|good or 2:blue=left|good
%        seqdir - sequence direction => 1:left or 2:right
%          resp - previous response => 1:left or 2:right (only if actor)
%
%  seqcat is equal to the sequence category (1:orange or 2:blue).
%
%  Valentin Wyart <valentin.wyart@ens.fr> - 09/2015

if nargin < 1
    error('Missing input arguments!');
end

if nargin < 3
    error('Invalid list of input arguments!');
end
seqcat = 1+(epimap ~= seqdir);


end