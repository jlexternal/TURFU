subject = 15; % choose from 1 through #
filename = dir(sprintf('Data/TURFU_S%02d_*.mat', subject));
filename = filename.name;
load(sprintf('Data/%s',filename));