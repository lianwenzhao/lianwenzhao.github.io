function [rmse, probHit2] = poissTest(poissWeight, featuresTest, trueCountTest, upBound)
U = length(featuresTest);
probHit2 = zeros(size(trueCountTest));
mse2 = zeros(1, U);
rmse = zeros(1,U);
for u=1:U
    poissMean = exp(poissWeight (1)+featuresTest{u}*poissWeight(2:end));
    if nargin>3
        poissMean = min(poissMean, upBound);
    end;
    probHit2(:, u) = exp(-poissMean);
    for i = 1:length(poissMean)
        probHit2(i, u) = probHit2(i, u) * poissMean(i)^trueCountTest(i,u) / factorial(trueCountTest(i,u));
    end;
    for i=1:length(poissMean)
        tmp = 0;
        count = 0;
        while (cdf('poiss', count, poissMean(i)) < 0.9999)
            tmp = tmp + pdf('poiss', count, poissMean(i)) * (count - trueCountTest(i,u))^2;
            count = count + 1;
        end;
        count = count + 1;
        tmp = tmp + pdf('poiss', count, poissMean(i)) * (count - trueCountTest(i,u))^2;
        mse2(u) = mse2(u) + tmp;
    end;
    rmse(u) = sqrt(mse2(u)/length(poissMean));
end;

