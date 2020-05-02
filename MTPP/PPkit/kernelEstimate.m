function [rate, tLocation, hOpt, mseAll] = kernelEstimate(events, Tstart, Tend, tLocation, deltaGrid, Tint, isTest, hCand)
%KERNELESTIMATE Summary of this function goes here
tLocation = unique([tLocation, events(events(:,1)>=Tstart,1)']); %%only works one mark
setDefault('hCand', deltaGrid*[2:2:((Tend-Tstart)/deltaGrid/10)]);
M = size(events,2);
if isTest
    backHist = max(floor(Tint/deltaGrid));
    tLocation = [Tstart-[backHist:-1:1]*deltaGrid , tLocation];
end;
rate = zeros(length(tLocation), M);
for m=1:M
    evt = events(:,m);
    stream.evt = evt;
    stream.Tstart = Tstart;
    stream.Tend = Tend;
    stream.N = length(stream.evt);
    stream.T = stream.Tend - stream.Tstart;
    if length(hCand)>1
        [hOpt, mse] = getKernelWidth(hCand, stream);
    else
        hOpt = hCand;
        mse = [];
    end;
    i = 0;
    for tt = tLocation
        i = i+1;
%         ts = max(Tstart, tt - hOpt);
%         te = min(Tend, tt + hOpt);
%         rate(i) = sum((evt>=(tt-hOpt))&(evt<=(tt+hOpt)));
        rate(i) = sum((evt>=(tt-hOpt))&(evt<=(tt)));%%%only consider events backward in time
    end;
    rate = rate/hOpt;%/2
    mseAll = [hCand; mse];
end
return;

function kh = secondMom(h, stream)
tmp = 0;
for i=1:stream.N
    idx =[1:1:stream.N];
    idx(i) =[];
    mid = abs(stream.evt(idx)- stream.evt(i));
    idxTmp = find(mid<=h);
    tmp = tmp + length(idxTmp);
    thr = min(stream.evt(i)-stream.Tstart, stream.Tend-stream.evt(i));
    tmp = tmp + sum(mid(idxTmp)>thr);
end;
kh = stream.T/stream.N^2*tmp;
return;

function [hOpt, mse] = getKernelWidth(hCand, stream)
numH = length(hCand);
kh = zeros(1, numH*2);
mu = stream.N/stream.T;
hCandDouble = [hCand, hCand + hCand(end)];
i = 0;
for h=hCandDouble
    i = i+1;
    kh(i) = secondMom(h, stream);
end;
dh = diff(hCand(1:2));
khIntegral = cumsum(([0,kh(1:end-1)]+kh)/2) * dh;
idx1 = [1:1:numH];
idx2 = [2:2:2*numH];
mse = (1-2*mu*kh(idx1))/2/mu./hCand + 1/4./hCand.^2.*khIntegral(idx2);
[~, hIdx] = min(mse);
hOpt = hCand(hIdx);



    
    
