%% param comparison
clc;

h   = [];
sig = [];
eta = [];
detRGB = [0, 0.4470, 0.7410];
infRGB = [0.8500, 0.3250, 0.0980];
selRGB = [0.4660, 0.6740, 0.1880];
load('selNoiseParamsLog');
post = selNoiseParamsLog;
load('selNoiseParamsLog_pred');
pred = selNoiseParamsLog;

%figure(1);
%hold on;
if false
for i = 1:30
    if isnan(post(i).past(1))
        continue;
    else
        h   = vertcat(h, [post(i).past(1) pred(i).future(1)]);
        sig = vertcat(sig, [post(i).past(2) pred(i).future(2)]);
        eta = vertcat(eta, [post(i).past(3) pred(i).future(3)]);
        %scatter(pred(i).future(1)-post(i).past(1), pred(i).future(3)-post(i).past(3),'.');
        %text(pred(i).future(1)-post(i).past(1), pred(i).future(3)-post(i).past(3), num2str(i));
    end
end
%plot([-1:1],[-1:1]);
%xlim([-.8 .8]);
%xlabel('H_{prediction} - H_{postdiction}');
%ylabel('\eta_{prediction} - \eta_{postdiction}');
%yline(0,':');
%xline(0,':');
%hold off;
figure(10);
violinplot([h], 'H','ViolinColor',detRGB,'EdgeColor',[1 1 1])
figure(11);
violinplot([sig], '\sigma','ViolinColor',infRGB,'EdgeColor',[1 1 1])
figure(12);
violinplot([eta], '\eta','ViolinColor',selRGB,'EdgeColor',[1 1 1])
end

% diff in sig to diff in H plot
postParams = [];
predParams = [];
X = [];
Y = [];
figure(21);
hold on;
for i = 1:30
    if isnan(post(i).past(2))
        continue;
    else
        postParams = vertcat(postParams, post(i).past(1:3));
        predParams = vertcat(predParams, pred(i).future(1:3));
        X = vertcat(X, pred(i).future(1)-post(i).past(1));
        Y = vertcat(Y, pred(i).future(2)-post(i).past(2));
        %text(pred(i).future(2)-post(i).past(2), pred(i).future(3)-post(i).past(3), num2str(i));
    end
end
scatter(X,Y,200,'.');
lsline;
[rsqrd pval1]  = corr(X,Y)
%[rho pval2] = corr(postParams, predParams)

plot([-1:1],[-1:1]);
yline(0,':');
xline(0,':');
xlim([-.5 .5]);
ylim([-.5 .5]);
xlabel('H_{prediction} - H_{postdiction}');
ylabel('\sigma_{prediction} - \sigma_{postdiction}');
hold off;

% diff in eta to diff in H plot
postParams = [];
predParams = [];
X = [];
Y = [];
figure(22);
hold on;
for i = 1:30
    if isnan(post(i).past(2))
        continue;
    else
        postParams = vertcat(postParams, post(i).past(1:3));
        predParams = vertcat(predParams, pred(i).future(1:3));
        X = vertcat(X, pred(i).future(1)-post(i).past(1));
        Y = vertcat(Y, pred(i).future(3)-post(i).past(3));
        %text(pred(i).future(2)-post(i).past(2), pred(i).future(3)-post(i).past(3), num2str(i));
    end
end
scatter(X,Y,200,'.');
lsline;
[rsqrd pval1]  = corr(X,Y)
%[rho pval2] = corr(postParams, predParams)

plot([-1:1],[-1:1]);
yline(0,':');
xline(0,':');
xlim([-.5 .5]);
ylim([-.5 .5]);
xlabel('H_{prediction} - H_{postdiction}');
ylabel('\eta_{prediction} - \eta_{postdiction}');
hold off;

% diff in sig to diff in eta plot
postParams = [];
predParams = [];
X = [];
Y = [];
figure(23);
hold on;
for i = 1:30
    if isnan(post(i).past(2))
        continue;
    else
        postParams = vertcat(postParams, post(i).past(1:3));
        predParams = vertcat(predParams, pred(i).future(1:3));
        X = vertcat(X, pred(i).future(2)-post(i).past(2));
        Y = vertcat(Y, pred(i).future(3)-post(i).past(3));
        %text(pred(i).future(2)-post(i).past(2), pred(i).future(3)-post(i).past(3), num2str(i));
    end
end
scatter(X,Y,200,'.');
lsline;
[rsqrd pval1]  = corr(X,Y)
%[rho pval2] = corr(postParams, predParams)
plot([-1:1],[-1:1]);
yline(0,':');
xline(0,':');
xlim([-.8 .8]);
ylim([-.8 .8]);
xlabel('\sigma_{prediction} - \sigma_{postdiction}');
ylabel('\eta_{prediction} - \eta_{postdiction}');
hold off;

if true
% compare params and plot
[hh hp hci hstats] = ttest(h(:,1), h(:,2))
[sigh sigp sigci sigstats] = ttest(sig(:,1), sig(:,2))
[etah etap etaci cetastats] = ttest(eta(:,1), eta(:,2))
figure(3);
hold on;
bar([1:2],[mean(h(:,1)) mean(h(:,2))],'FaceColor',detRGB);
errorbar([mean(h(:,1)) mean(h(:,2))],[std(h(:,1)) std(h(:,2))],'LineStyle','none','LineWidth',2,'CapSize',0,'Color','k')
xticks([1:2]);
xticklabels({'postdiction' 'prediction'});
ylabel('value of H');
hold off;
figure(4);
hold on;
bar([1:2],[mean(sig(:,1)) mean(sig(:,2))],'FaceColor',infRGB);
errorbar([mean(sig(:,1)) mean(sig(:,2))], [std(sig(:,1)) std(sig(:,2))],'LineStyle','none','LineWidth',2,'CapSize',0,'Color','k')
xticks([1:2]);
xticklabels({'postdiction' 'prediction'});
ylabel('value of \sigma');
hold off;
figure(5);
hold on;
bar([1:2],[mean(eta(:,1)) mean(eta(:,2))],'FaceColor',selRGB);
errorbar([mean(eta(:,1)) mean(eta(:,2))], [std(eta(:,1)) std(eta(:,2))],'LineStyle','none','LineWidth',2,'CapSize',0,'Color','k');
xticks([1:2]);
xticklabels({'postdiction' 'prediction'});
ylabel('value of \eta');
hold off;
end

%% alternative confidence thing

n_subjects = 30;
excluded = [14 20 21 22 27];

avgConf = zeros(25,1);
postConfs = zeros(25,8);
predConfs = zeros(25,8);
subjCtr = 1;
for isubj = 1:n_subjects
    
    if ismember(isubj, excluded)
       continue; 
    end
    
    subject  = isubj;
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    % total cross-condition average
    avgConf(subjCtr) = 0;
    for iblock = 3:10
        hiConf = expe.rslt(iblock).conf(2:73);
        hiConf = hiConf(hiConf==2);
        avgConf(subjCtr) = avgConf(subjCtr) + numel(hiConf);
    end
    nTotalTrials = 72*8;
    avgConf(subjCtr) = avgConf(subjCtr)./nTotalTrials;
    
    postNRev = 0;
    predNRev = 0;
    
    for iblock = 3:10
    nRevs = 0;
    tempConf = zeros(1,8);
        % loop through the trials 
        for i = 3:72    % start from 3 since starting from 2 will compare with the null choice
            % REVERSAL CURVE analysis
            % identify trials after switch 
            if expe.blck(iblock).seqdir(i) ~= expe.blck(iblock).seqdir(i-1)
                nRevs = nRevs + 1;
                % get subject response and confidence
                inRevs = -3;
                while inRevs ~= 5
                    if expe.rslt(iblock).conf(i+inRevs) == 2
                        tempConf(inRevs+4) = tempConf(inRevs+4) + 1;
                    end
                    inRevs = inRevs + 1;
                end
            end
        end
        
        %store by condition
        if expe.blck(iblock).taskid == 1
            postNRev = postNRev + nRevs;
            postConfs(subjCtr,:) = postConfs(subjCtr,:) + tempConf;
        else
            predNRev = predNRev + nRevs;
            predConfs(subjCtr,:) = predConfs(subjCtr,:) + tempConf;
        end
    end
    postConfs(subjCtr,:) = postConfs(subjCtr,:)./postNRev;
    predConfs(subjCtr,:) = predConfs(subjCtr,:)./predNRev;
    subjCtr = subjCtr + 1;
end

postConfs = postConfs - avgConf;    
predConfs = predConfs - avgConf;

%% plotting
postRGB = [1 .75 .5];
predRGB = [.75 .75 1];

figure(1);
hold on;
errorbar([1:8], mean(postConfs),std(postConfs)/sqrt(length(postConfs)),'LineStyle','none',...
                    'Color',postRGB);    % condition 1 stdev over single subject means
scatter([1:8], mean(postConfs),100,'filled','MarkerFaceColor',postRGB);
errorbar([1:8], mean(predConfs),std(predConfs)/sqrt(length(predConfs)),'LineStyle','none',...
                    'Color',predRGB);    % condition 1 stdev over single subject means
scatter([1:8], mean(predConfs),100,'filled','d','MarkerFaceColor',predRGB);
ylim([-.2 .1]);
yline(0);

hold off;



%% script to make matrix for Valentin

dataMat = struct;

n_subjects = 30;
excluded = [14 20 21 22 27];
subjCtr = 1;
for isubj = 1:n_subjects
    
    if ismember(isubj, excluded)
       continue; 
    end
    
    dataMat(subjCtr).blck = struct;
    
    %load subject information
    subject  = isubj;
    filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
    filename = filename.name;
    load(sprintf('Data/%s',filename));
    
    blocks1 = find([expe.blck.taskid] == 1 & [expe.blck.condtn]~=3);
    blckCtr = 1;
    for block = blocks1
        for itrial = 2:73
            if expe.rslt(block).resp(itrial) == expe.rslt(block).resp(itrial-1)
                dataMat(subjCtr).blck(blckCtr).stay(itrial-1) = 1;
                dataMat(subjCtr).blck(blckCtr).signedllr(itrial-1) = abs(expe.blck(block).seqllr(itrial-1));
            else
                dataMat(subjCtr).blck(blckCtr).stay(itrial-1) = 0;
                dataMat(subjCtr).blck(blckCtr).signedllr(itrial-1) = -abs(expe.blck(block).seqllr(itrial-1));
            end
            
        end
        blckCtr = blckCtr + 1;
    end
    
    blocks2 = find([expe.blck.taskid] == 2 & [expe.blck.condtn]~=3);
    for block = blocks2
        for itrial = 2:73
            if expe.rslt(block).resp(itrial) == expe.rslt(block).resp(itrial-1)
                dataMat(subjCtr).blck(blckCtr).stay(itrial-1) = 1;
                dataMat(subjCtr).blck(blckCtr).signedllr(itrial-1) = abs(expe.blck(block).seqllr(itrial-1));
            else
                dataMat(subjCtr).blck(blckCtr).stay(itrial-1) = 0;
                dataMat(subjCtr).blck(blckCtr).signedllr(itrial-1) = -abs(expe.blck(block).seqllr(itrial-1));
            end
            
        end
        blckCtr = blckCtr + 1;
    end
    subjCtr = subjCtr + 1;
end