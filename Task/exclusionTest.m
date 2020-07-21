% Testing subject responses against random answers on the past condition
%
%   Name:   exclusionTest.m
% Author:   Jun Seok Lee
%   Date:   March 2019
%   Type:   Script
%   Note:   For TURFU experiment

n_subjects  = 30;
excluded    = [21]%[8 14 21]; % Subject numbers to exclude

losers      = zeros(10,30);
n_runs      = 100;


%% Compare subject responses to random answers
for runs = 1:n_runs
    for isubj = 1:n_subjects

        if ismember(isubj,excluded)
            continue;
        end

        % load subject data
        filename = dir(sprintf('Data/TURFU_S%02d_*.mat', isubj));
        filename = filename.name;
        load(sprintf('Data/%s',filename));

        % go through real blocks and compare
        for i = 3:10

            %if expe.blck(i).taskid == 2 % do not compare results from future condition
             %   continue;
            %end
            
            randWrong   = expe.blck(i).seqdir(1:72) - randi(2,1,72);
            randWrong   = numel(randWrong(randWrong~=0));

            subjWrong   = expe.rslt(i).resp(2:73) - expe.blck(i).seqdir(1:72);
            subjWrong   = numel(subjWrong(subjWrong~=0));
            if subjWrong > randWrong
                disp(['Subject ' num2str(isubj) ' performed worse than random on block ' num2str(i)]);
                losers(i, isubj) = losers(i,isubj)+1;
            end
        end
    end
end

%% List of future conditions
for isubj = 1:n_subjects

    if ismember(isubj,excluded)
        continue;
    end

    % load subject data
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', isubj));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    disp('');
    disp(['Subject ' num2str(isubj)]);
    
    for i = 3:10
        if expe.blck(i).taskid == 2
            disp(i);
        end
    end

end
%% Plots
for i = 3:10
    subplot(8,1,i-2);
    bar([1:30], losers(i,:)/n_runs);
    xticks([1:30]);
    ylim([0 1]);
    title(sprintf('Percentage of worse than random performance after %i runs in Block %i ', n_runs,i-2));
    hold on;
    plot(excluded, zeros(1,numel(excluded)),'x');
    hold off;
end