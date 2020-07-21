%% Load and analyzed data for TURFU
% comput both P_stay and reversal curves
clear all
%clc
%close all
subjlist = setdiff(1:30,[14 20 21 22 27]);
%subjlist =setdiff(1:30,[14,21]);

nsubj = numel(subjlist);

datapath = '../Data';
addpath('~/Documents/aurelien/actobs/script2')

fparam = 'woe'; % fixed parameter for both condition: lps (lapse is fixed), pse (shift is fixed) or woe (slope is fixed)
col = @(x)x(:);

for isubj = 1:nsubj
    subj = subjlist(isubj);
    
    fprintf('processing subj %02d\n',subj)
    % load raw behavior
    str_4cmp = sprintf('S%02d',subj);
    f = dir(datapath);
    for i = 4:size(f,1)
        if strcmp(f(i).name(7:9),str_4cmp)
            break
        end
    end
    
    dataname = sprintf('%s/%s',datapath,f(i).name);
    load(dataname)
    % ignore training blocks
    blck = expe.blck(3:10);
    rslt = expe.rslt(3:10);
    
    % extract relevant information relative to task
    taskid_resp = kron([blck.taskid],ones(1,72));
    taskid_resp = taskid_resp(:);
    
    epinum = [];
    epipos = [];
    seqllr = [];
    seqdir = [];
    
    for irun = 1:8
        seqdir = [seqdir,blck(irun).seqdir];
        epipos = [epipos,blck(irun).epipos]; % episode position within a reversal
        seqllr = [seqllr,blck(irun).seqllr];
        epinum = [epinum,blck(irun).epinum]; % episode position within a run
    end
    
    seqind = repmat(1:72,1,8);
    % extract relevant information relative to behavior
    resp   = [];
    rd_bef = [];
    rt     = [];
    
    for irun = 1:8
        resp   = [resp, rslt(irun).resp(2:end)];
        rd_bef = [rd_bef, rslt(irun).resp(1:end-1)];
        rt     = [rt, rslt(irun).rt(2:end)];
    end
        
        ifilt  = seqind > 1;
        % behavioral analyses: 1. P(stay) = f(cnf)
        xcnf = col(seqllr(ifilt).*(3-2*rd_bef(ifilt)));
        y    = col(resp(ifilt) == rd_bef(ifilt));
        xvec = (-4:0.01:+4)';
        
        taskid_resp = taskid_resp(seqind > 1);
        switch fparam
            case 'pse' % fixed pse
                [pval,lps,llh] = logreg_psewoe_condition_pse(xcnf,y,'logit',taskid_resp);
                
                for itask = 1:2
                % compute fit
                pse_task_resp(isubj,itask) = fzero(@(p)lps(itask)+(1-lps(itask))./(1+exp(-[1,p]*pval(:,itask)))-0.5,[-10,+10]);
                presp_resp(isubj,:,itask) = lps(itask)+(1-lps(itask))./(1+exp(-cat(2,ones(size(xvec)),xvec)*pval(:,itask)));               
                woe_all(isubj,itask) = pval(2,itask);
                lps_all(isubj,:) = lps;
                end
                llh_all(isubj) = llh;
                
             %  aic_subj(isubj,itask) = aic;
              %  bic_subj(isubj,itask) = bic;
            case 'woe' % fixed woe
                 [pval,lps,llh] = logreg_psewoe_condition_woe(xcnf,y,'logit',taskid_resp);
                
                for itask = 1:2
                % compute fit
                pse_task_resp(isubj,itask) = fzero(@(p)lps(itask)+(1-lps(itask))./(1+exp(-[1,p]*pval(:,itask)))-0.5,[-10,+10]);
                presp_resp(isubj,:,itask) = lps(itask)+(1-lps(itask))./(1+exp(-cat(2,ones(size(xvec)),xvec)*pval(:,itask)));               
                woe_all(isubj,itask) = pval(2,itask);
                end
                lps_all(isubj,:) = lps;
                llh_all(isubj) = llh;
                
                     case 'lps' % fixed pse
                [pval,lps,llh] = logreg_psewoe_condition_lps(xcnf,y,'logit',taskid_resp);
                
                for itask = 1:2
                % compute fit
                pse_task_resp(isubj,itask) = fzero(@(p)lps+(1-lps)./(1+exp(-[1,p]*pval(:,itask)))-0.5,[-10,+10]);
                presp_resp(isubj,:,itask) = lps+(1-lps)./(1+exp(-cat(2,ones(size(xvec)),xvec)*pval(:,itask)));               
                woe_all(isubj,itask) = pval(2,itask);
                
                lps_all(isubj,itask) = lps;
                end
                 llh_all(isubj) = llh;
                
             %  aic_subj(isubj,itask) = aic;
              %  bic_subj(isubj,itask) = bic;
        end
                
        for itask = 1:2
        ifilt = taskid_resp == itask;
        % behavioral analyses: 1. P(stay) = f(cnf)
        xcnf = col(seqllr(ifilt).*(3-2*rd_bef(ifilt)));
        y    = col(resp(ifilt) == rd_bef(ifilt));
        
                bin_type = 'fixed';
        switch bin_type
            case 'fixed'
               xmin = [-inf,  -3,  -2,  -1,   0,  +1,  +2,  +3];
               xmax = [  -3,  -2,  -1,   0,  +1,  +2,  +3,+inf];
                
               %  xmin = [-inf, -4 ,-3,  -2,  -1,   0,  +1,  +2,  +3 +4];
               %  xmax = [  -4,-3,  -2,  -1,   0,  +1,  +2,  +3, +4,+inf];
                nbin = numel(xmax);
                for i = 1:nbin
                    respfilt = xcnf > xmin(i) & xcnf < xmax(i);
                    xbin_resp(isubj,i,itask) = -4.5+i;
                    ybin_resp(isubj,i,itask) = mean(y(respfilt));
                    
                end
            otherwise
                error('undefined binning type!');
        end
        end        

end

%% Plot P(stay) = f(evidence conf)
yavg = squeeze(mean(presp_resp,1));
yerr = squeeze(std(presp_resp,[],1)./sqrt(nsubj));
ybin = squeeze(mean(ybin_resp,1));
ebin = squeeze(std(ybin_resp,[],1)./sqrt(nsubj));

figure
hold on 
%rgb = cat(1,[0.5,0.75,1],[1,0.5,0.75]);
rgb = cat(1,[1,0.75,0.5],[0.75,0.75,1])
% fitted curves
patch([xvec',fliplr(xvec')],[yavg(:,1)'+yerr(:,1)',fliplr(yavg(:,1)'-yerr(:,1)')],0.25*rgb(1,:)+0.75,'EdgeColor','none')
patch([xvec',fliplr(xvec')],[yavg(:,2)'+yerr(:,2)',fliplr(yavg(:,2)'-yerr(:,2)')],0.25*rgb(2,:)+0.75,'EdgeColor','none')
plot(xvec,yavg(:,1),'-','LineWidth',2,'Color',rgb(1,:))
plot(xvec,yavg(:,2),'-','LineWidth',2,'Color',rgb(2,:))

% avg bin
plot(xbin_resp(1,:,1),ybin(:,1),'o','MarkerEdgeColor',rgb(1,:),'MarkerFaceColor',rgb(1,:))
plot(xbin_resp(1,:,1),ybin(:,2),'o','MarkerEdgeColor',rgb(2,:),'MarkerFaceColor',rgb(2,:))

for i = 1:length(xbin_resp(1,:,1))
    plot([xbin_resp(1,i,1),xbin_resp(1,i,1)], [ybin(i,1) + ebin(i,1),ybin(i,1) - ebin(i,1)],'Color',rgb(1,:))
    plot([xbin_resp(1,i,1),xbin_resp(1,i,1)], [ybin(i,2) + ebin(i,2),ybin(i,2) - ebin(i,2)],'Color',rgb(2,:))
end

pbar = 1.5;
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
%et(gca,'Xtick',[0 0.2 0.4 0.6],'XtickLabel',{'0','0.2','0.4','0.6'},'Xlim',[-0.1,0.7])
set(gca,'Ytick',0:0.2:1,'YtickLabel',{'0','0.2','0.4','0.6','0.8','1'})
xlabel('evidence for stay (log-odds)','FontSize',8)
ylabel('p(stay)','FontSize',8)
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
%title(pop,'FontName','Helvetica','FontSize',7.2)

%% Barplot of parameters and stats : WOE
[~,p,~,t] = ttest(woe_all(:,2),woe_all(:,1))
rgb = cat(1,[1,0.75,0.5],[0.75,0.75,1]);
xavg = squeeze(mean(woe_all,1));
xerr = squeeze(std(woe_all,[],1)./sqrt(nsubj));

figure;
hold on
for itask = 1:2
    bar(itask,xavg(itask),'FaceColor',rgb(itask,:),'LineWidth',1.5)
    plot([itask,itask],[xavg(itask) + xerr(itask), xavg(itask) - xerr(itask)],'k-','LineWidth',2)
    
end
xlim([0,3])
xticks([1,2])
xticklabels({'Cue','Outcome'})
ylim([0,2])
pbar = 3/4;
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
%et(gca,'Xtick',[0 0.2 0.4 0.6],'XtickLabel',{'0','0.2','0.4','0.6'},'Xlim',[-0.1,0.7])
set(gca,'Ytick',0:1:2,'YtickLabel',{'0','1','2'})
%xlabel('evidence for stay (log-odds)','FontSize',8)
ylabel('weight of evidence (a.u)','FontSize',8)
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
%title(pop,'FontName','Helvetica','FontSize',7.2)

%% Barplot of parameters and stats : PSE
[~,p,~,t] = ttest(pse_task_resp(:,2),pse_task_resp(:,1))
rgb = cat(1,[1,0.75,0.5],[0.75,0.75,1]);
xavg = squeeze(mean(-pse_task_resp,1));
xerr = squeeze(std(pse_task_resp,[],1)./sqrt(nsubj));

figure;
hold on
for itask = 1:2
    bar(itask,xavg(itask),'FaceColor',rgb(itask,:),'LineWidth',1.5)
    plot([itask,itask],[xavg(itask) + xerr(itask), xavg(itask) - xerr(itask)],'k-','LineWidth',2)
    
end
xlim([0,3])
xticks([1,2.25])
xticklabels({'Cue','Outcome'})
ylim([0,3])
pbar = 3/4;
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
%et(gca,'Xtick',[0 0.2 0.4 0.6],'XtickLabel',{'0','0.2','0.4','0.6'},'Xlim',[-0.1,0.7])
set(gca,'Ytick',[0,0.75,1.5,2.25],'YtickLabel',{'0','0.75','1.5','2.25'})
%xlabel('evidence for stay (log-odds)','FontSize',8)
ylabel('point of subjective equivalence (a.u)','FontSize',8)
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
%title(pop,'FontName','Helvetica','FontSize',7.2)

%% proportion stay
[~,p,~,t] = ttest(lps_all(:,2),lps_all(:,1))
rgb = cat(1,[1,0.75,0.5],[0.75,0.75,1]);
xavg = squeeze(mean(lps_all,1));
xerr = squeeze(std(lps_all,[],1)./sqrt(nsubj));

figure;
hold on
for itask = 1:2
    bar(itask,xavg(itask),'FaceColor',rgb(itask,:),'LineWidth',1.5)
    plot([itask,itask],[xavg(itask) + xerr(itask), xavg(itask) - xerr(itask)],'k-','LineWidth',2)
    
end
xlim([0,3])
xticks([1,2])
xticklabels({'Cue','Outcome'})
ylim([0,0.15])
pbar = 3/4;
set(gca,'Layer','top','Box','off','PlotBoxAspectRatio',[pbar,1,1]);
%et(gca,'Xtick',[0 0.2 0.4 0.6],'XtickLabel',{'0','0.2','0.4','0.6'},'Xlim',[-0.1,0.7])
set(gca,'Ytick',[0,0.1],'YtickLabel',{'0','0.1'})
%xlabel('evidence for stay (log-odds)','FontSize',8)
ylabel('negative lapse (a.u)','FontSize',8)
set(gca,'TickDir','out','TickLength',[1,1]*0.02/max(pbar,1));
set(gca,'FontName','Helvetica','FontSize',7.2);
%title(pop,'FontName','Helvetica','FontSize',7.2)

%clc
mean(llh_all)
% mean llh fixed pse = -210.1847
% mean llh fixed lps = -210.0820 wins,
% mean llh fixed woe = -212.0716

%% mixed anova: evidence x condition
[~,p,~,t] = ttest(ybin_resp(:,:,2),ybin_resp(:,:,1))
dm = reshape(ybin_resp,[size(ybin_resp,1).*size(ybin_resp,2),size(ybin_resp,3)]);

datamat = cat(2,dm(:,2),dm(:,1));
databtw = [];
for i = 1:8
databtw = cat(1,databtw,i*ones(nsubj,1));
end
[tbl,rm] = simple_mixed_anova(datamat,databtw,{'task'},{'evidence'});
