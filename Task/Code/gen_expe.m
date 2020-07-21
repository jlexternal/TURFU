function [expe] = gen_expe(suj,taskeq)
%  GEN_EXPE  Generate ACTOBS experiment
%
%  Usage: [expe] = GEN_EXPE(suj,taskeq)
%
%  where suj is the subject number. Experiment parameters are counter-balanced
%  after groups of 8 subjects. The experiment structure expe contains the block
%  substructure blck required to run the experiment.
%
%  The optional parameter taskeq controls the degree of equalization between
%  the two tasks (true for full theoretical equalization between the two tasks,
%  false for partial equalization between them). The default is true.
%
%  Valentin Wyart <valentin.wyart@ens.fr> - 09/2015

% check input arguments
if nargin < 2
    % apply full equalization between tasks
    taskeq = true;
end
if nargin < 1
    error('Missing subject number!');
end
% task identifier => 1:observer or 2:actor
% 1st bit of subject number
if bitget(suj,1)
    taskid = [1 2  1 1 2 2  1 1 2 2];
else
    taskid = [2 1  2 2 1 1  2 2 1 1];
end

% task condition => 1:stable or 2:volatile or 3:practice
% 2nd bit of subject number
if bitget(suj,2)
    condtn = [3 3  1 2 1 2  2 1 2 1];
else
    condtn = [3 3  2 1 2 1  1 2 1 2];
end

% episode mapping => 1:pink=good|left or 2:blue=good|left
% 3rd bit of subject number
if bitget(suj,3)
    epimap = [1 1  2 1 2 1  2 1 2 1];
else
    epimap = [2 2  1 2 1 2  1 2 1 2];
end
nblck = length(taskid);

% initial episode direction => 1:left or 2:right
% pseudo-randomized independently for each subject
epidir = kron([1,2],ones(1,nblck/2));
while true
    epidir = epidir(randperm(nblck));
    if ~HasConsecutiveValues(epidir,3)
        break
    end
end

if taskeq
    % full equalization between tasks
    isub = find(taskid == taskid(1)); % block indices of 1st task
    jsub = find(taskid ~= taskid(1)); % block indices of 2nd task
    % check ordering of conditions
    if any(condtn(isub) ~= condtn(jsub) | epimap(isub) ~= epimap(jsub))
        error('Invalid ordering of conditions!');
    end
    for i = 1:length(isub)
        % generate block for 1st task
        iblck = isub(i);
        cfg        = [];
        cfg.condtn = condtn(iblck);
        cfg.epimap = epimap(iblck);
        cfg.epidir = epidir(iblck);
        blck = gen_blck(cfg);
        blck.taskid = taskid(iblck);
        expe.blck(iblck) = blck;
        % copy block information for 2nd task
        jblck = jsub(i);
        expe.blck(jblck) = expe.blck(iblck);
        expe.blck(jblck).taskid = taskid(jblck);
        if epidir(jblck) ~= epidir(iblck)
            % flip sequence directions and evidence
            expe.blck(jblck).seqdir = 3-expe.blck(jblck).seqdir;
            expe.blck(jblck).seqllr = -expe.blck(jblck).seqllr;
        end
    end
else
    % partial equalization between tasks
    for iblck = 1:nblck
        % generate block
        cfg        = [];
        cfg.condtn = condtn(iblck);
        cfg.epimap = epimap(iblck);
        cfg.epidir = epidir(iblck);
        blck = gen_blck(cfg);
        blck.taskid = taskid(iblck);
        expe.blck(iblck) = blck;
    end
end

% set taskid as 1st field in block substructure
n = length(fieldnames(expe.blck));
expe.blck = orderfields(expe.blck,[n,(1:n-1)]);

end