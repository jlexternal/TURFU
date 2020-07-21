function [repcat] = get_repcat(taskid,epimap,resp)
%  GET_REPCAT  Get response category
%
%  Usage: [repcat] = GET_REPCAT(taskid,epimap,[resp])
%
%  where taskid - task identifier => 1:observer or 2:actor
%        epimap - episode mapping => 1:orange=left|good or 2:blue=left|good
%          resp - response => 1:left or 2:right (only if observer task)
%
%  repcat is equal to the chosen category in the observer condition (1:orange
%  or 2:blue), and equal to the requested category in the actor condition (and
%  is therefore independent from the laterality of the response).
%
%  Valentin Wyart <valentin.wyart@ens.fr> - 09/2015

if nargin < 1
    error('Missing input arguments!');
end

if taskid == 1 % observer task
    if nargin < 3
        error('Invalid list of input arguments!');
    end
    repcat = 1+(epimap ~= resp);
elseif taskid == 2 % actor task
    if nargin < 2
        error('Invalid list of input arguments!');
    end
    % subject is by definition trying to generate the requested category
    % repcat is therefore independent from resp
    repcat = epimap;
else
    error('Invalid task identifier!');
end

end