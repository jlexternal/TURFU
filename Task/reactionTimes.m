excluded = [14 20 21 22 27];
n_subjects = 30;
n_trialsTot = 72*4;
n_conditions = 2;
rts = zeros(n_trialsTot,n_subjects-numel(excluded),n_conditions); % (number of subjects, 288 trials, 2 conditions)

isubjCtr = 1;
for isubj = 1:n_subjects
    
    if ismember(isubj, excluded)
       continue; 
    end
    
    %load subject information
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', isubj));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    itrialCtr1 = 1;
    itrialCtr2 = 1;
    for iblck = 3:10
        if expe.blck(iblck).taskid == 1
            rts(itrialCtr1:itrialCtr1+71,isubjCtr,1)= expe.rslt(iblck).rt(2:73).';
            itrialCtr1 = itrialCtr1 + 72;
        elseif expe.blck(iblck).taskid == 2
            rts(itrialCtr2:itrialCtr2+71,isubjCtr,2)= expe.rslt(iblck).rt(2:73).';
            itrialCtr2 = itrialCtr2 + 72;
        end
    end
    isubjCtr = isubjCtr + 1;
end

rtmeans1 = mean(rts(:,:,1));
rtmeans2 = mean(rts(:,:,2));

mean(rtmeans1)
mean(rtmeans2)

ttest(rtmeans1,rtmeans2)



