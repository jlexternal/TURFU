% script to make paramLog structure not as stupid as I first coded it
% Author: Jun Seok Lee

load('paramLog','paramLog');

paramsLog   = struct;
n_subjs     = numel(paramLog.dir.past(:,1));
for i = 1:n_subjs
    paramsLog(i).past   = paramLog.dir.past(i,:);
    %paramsLog(i).future = paramLog.dir.future(i,:);
    %paramsLog(i).low    = paramLog.vol.low(i,:);
    %paramsLog(i).high   = paramLog.vol.high(i,:);
end

%clearvars -except paramsLog

save('paramsLog_test');

%% Fix glazeOrigParamLog
load('glazeOrigParamLog','glazeOrigParamLog');

glazeOrigParamsLog   = struct;
n_subjs     = numel(glazeOrigParamLog.dir.past);
for i = 1:n_subjs
    disp(num2str(i));
    glazeOrigParamsLog(i).past   = glazeOrigParamLog.dir.past(i);
    glazeOrigParamsLog(i).future = glazeOrigParamLog.dir.future(i);
    %glazeOrigParamsLog(i).low    = glazeOrigParamLog.vol.low(i);
    %glazeOrigParamsLog(i).high   = glazeOrigParamLog.vol.high(i);
end

clearvars -except glazeOrigParamsLog

save('glazeOrigParamsLog');