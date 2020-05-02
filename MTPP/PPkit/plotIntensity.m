function [ fig] = plotIntensity( trainTime, trainMean, trainStd, testTime, testMean, testStd )
%PLOTINTENSITY Summary of this function goes here
%% plot rate
fig = figure;
if nargin<4
    shadedErrorBar([trainTime], [trainMean], [trainStd],'k');
else
    shadedErrorBar([trainTime; testTime], [trainMean; testMean], [trainStd; testStd],'k');
end;
