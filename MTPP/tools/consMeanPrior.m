function g = consMeanPrior(markedEvents, Tstart, Tend)
U = length(markedEvents);
numEvent = 0;
for u=1:U
    numEvent = numEvent + size(markedEvents{u},1);
end
g = sqrt(numEvent/(Tend-Tstart)/U);

