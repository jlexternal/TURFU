%  RUNNER  Runner script for ACTOBS experiment
%
%  This script runs the ACTOBS experiment and saves automatically the obtained
%  results in the ../Data folder. This script should be run from its containing
%  folder (which contains all the functions of the ACTOBS experiment).
%
%  Note that the data will not be saved after exiting if an error was caught
%  during the experiment! Check the command window for errors.
%
%  Valentin Wyart <valentin.wyart@ens.fr> - 09/2015

% clear workspace
clear all java
close all hidden
clc

addpath('./Toolboxes/Rand');

% setup random number generator
SetupRand;

% get participant information
argindlg = inputdlg({'Subject number (two-digit)','Use eye-tracker setup? (yes/no)'},'ACTOBS',1,{'','no'});
if isempty(argindlg)
    error('experiment cancelled!');
end
hdr      = [];
hdr.suj  = str2num(argindlg{1});
hdr.date = datestr(now,'yyyymmdd-HHMM');

% use CENIR MEG setup?
eyetrack = strcmp(lower(argindlg{2}),'yes');

% generate experiment
expe = gen_expe(hdr.suj,true);

% add header to experiment structure
expe.hdr = hdr;
expe = orderfields(expe,{'hdr','blck'});
% run experiment
[expe,aborted,errmsg] = run_expe(expe,eyetrack);
if ~isempty(errmsg)
    rethrow(errmsg);
end

% save data
datapath = ['..',filesep,'Data'];
filename = sprintf('TURFU_S%02d_%s',hdr.suj,hdr.date);
if aborted
    filename = [filename,'_aborted'];
end
filename = [filename,'.mat'];
filename = fullfile(datapath,filename);
save(filename,'expe');
