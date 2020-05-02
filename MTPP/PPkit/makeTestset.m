function [ testSet, testTime ] = makeTestset( markedEvents, Tint, Tend, M )
mid = stream2instance(markedEvents, Tint, Tend, M);
U = length(mid);
testSet = cell(1,U);
testTime = cell(1, U);
for u= 1:U
    testSet{u} = mid{u}.feature;
    testTime{u} = [mid{u}.time(1)-mid{u}.duration(1); mid{u}.time];
end

