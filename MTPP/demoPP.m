%%%This software implements the algorithm in 
%%% Lian, Wenzhao, et al. "A multitask point process predictive model." 
%%% Proceedings of the 32nd International Conference on Machine Learning (ICML-15). 2015.
%%% Please contact wl89@duke.edu for questions.
clear; close all;
addpath(genpath('./'));
load('./PPkit/gz500.mat')
dataSource = 1; %%1 for toy data; 2 for real data
if (dataSource==1)
    %% generate toy data
    U = 10;%%number of users
    C = 1;%%number of clusters of users
    M = 1;%%number of marks
    Tstart= 0;
    Tend= 50;%100
    TtestEnd = 100;%200
    Tint= [5 1]';
    rateOption = 'pseudo';%'linear';
    ratePar = rateGenerator(rateOption, U);
    plotOption = 'on';
    model0.ratePar = ratePar;
    e0=cell(1,U);for u=1:U;e0{u}=[0 1];end;
    [ markedEventsAll, eventCountAll, trueInstAll, rateTrueAll] = forwardSampling(Tstart, TtestEnd, Tint, M, model0, e0, 1, plotOption);
    for u=1:U
        markedEventsTrain{u}= markedEventsAll{u}(markedEventsAll{u}(:,1)<=Tend,:);
        rateTrueTrain{u} = rateTrueAll{u}(rateTrueAll{u}(:,1)<=Tend,:);
        markedEventsTest{u}= markedEventsAll{u}(markedEventsAll{u}(:,1)>Tend,:);
        rateTrueTest{u} = rateTrueAll{u}(rateTrueAll{u}(:,1)>Tend,:);
    end;
    disp('data generated....');
elseif (dataSource==2)
    %% real data
    %%markedEventsAll=getRealData()...
    %%Parameter setting...
end;
%%%% Prepare data
[ pieceAll ] = stream2instance(markedEventsTrain, Tint, Tend, M);%%process data into feature instances
[ piecePop, pieceAll, piecePopRep ] = buildTable( pieceAll );
[N, D] = size(piecePop);
%% Inference
% Fix seeds
rand('seed', 2e5);
randn('seed', 2e5);

% Parameter setting.
[~, options.pseudo]=kmeans(piecePop,20);
options.m = size(options.pseudo,1);
options.prior.xi = 0.1;
options.prior.nu = options.m+2;
options.ratio = 10;
options.prior.g =  sqrt(options.ratio) * consMeanPrior(markedEventsTrain, Tstart, Tend);

% Inference mode setup
trops(1) = 1000; % number of iterations
trops(2) = 1; % type of optimizations
trops(3) = 0; %%optimize GP hyper only, 3; 0 opt nothing
%%% opt GP hyper and inducing variables, 1; opt inducing variables only, 2
trops(4) = 20; %number of EM iterations

% type of inducing variables
options.Likelihood = 'PointProcess';
options.indType = 'pseudoIns'; % 'pseudoIns' or 'weights'
options.objectFunc = 'var';
options.vg = 'on';

rateTotal = [];
deltaGrid=(Tend-Tstart)/200;
tLocation = [Tstart:deltaGrid:Tend];
hCand = deltaGrid*[1:1:20];
for u=1:U
    [rateKernel, tLocation, hOpt, mseH] = kernelEstimate(markedEventsTrain{u}(:,1), Tstart, Tend,  tLocation, deltaGrid, Tint, 0, hCand);
    rateTotal = [rateTotal; rateKernel];
end;
options.stdRate = std(rateTotal);
model = varsgpCreate('sum', piecePop, [], options, pieceAll);
Xuinit = model.Xuinit;

% initialization of the model hyperparameters
logtheta0(1:D+2,1) = model.GP.logtheta;
model.gz = gz;
model.dgz = dgz;

% train the model
[model, margLogL] = varsgpTrain(model, trops);

disp('Training finished...')

figure;plot(margLogL); title('Lower bound','FontSize',16);
set(gca,'FontSize',16);

%% prediction
TpredEnd = Tend + 5; 
MCiter = 100;
[ predEventAll, predCount, predInstAll, predRateAll] = forwardSampling(Tend, TpredEnd, Tint, M, model, markedEventsTrain, MCiter, 'off');
for u=1:U;
    trueCount(u) = sum((markedEventsAll{u}(:,1)>=Tend).*(markedEventsAll{u}(:,1)<TpredEnd));
end;
if dataSource==1
    model0.ratePar = ratePar;
    [ predEventAll0, predCount0, predInstAll0, predRateAll0] = forwardSampling(Tend, TpredEnd, Tint, M, model0, markedEventsTrain, MCiter, 'off');
    figure;
    Umax=5;
    for u=1:Umax
        subplot(Umax,1,u);
        [nn, pp] = hist([predCount(:,u) predCount0(:,u)]);
        bar(pp, nn/MCiter);
        hold on;plot([trueCount(u), trueCount(u)], ylim,'r','lineWidth',2);
        legend('Inferred', 'Truth','Observed count');
        KL(u)= sum(nn(:,2)/MCiter.*log(eps+(nn(:,2)/MCiter)./(eps+nn(:,1)/MCiter)));
    end;
end;
%% forward test for count in intervals
maxNumTest = 50;
MCiter = 100;
Tduration = 5;
[ rmse, probHit, featuresTest, trueCountTest ] = forwardTest( Tend, Tduration, TtestEnd, Tint, M, model, markedEventsAll, maxNumTest, MCiter);
if dataSource==1
    [ rmse0, probHit0 ] = forwardTest( Tend, Tduration, TtestEnd, Tint, M, model0, markedEventsAll, maxNumTest, MCiter);
end;

% poisson regression
Tstep = 0.5*Tduration;
[xInterval, yInterval] =  buildIntervalSet(markedEventsTrain, M, Tstart, Tend, Tint, Tduration,  Tstep);
xPoiss = cell2mat(xInterval');
yPoiss = cell2mat(yInterval');
poissWeight = glmfit(xPoiss, yPoiss, 'poisson');
upBound = max(yPoiss)*3;
[rmse2, probHit2] = poissTest(poissWeight, featuresTest, trueCountTest, upBound);

% display results
if dataSource==1;disp(['RMSE with TRUE: ' num2str(mean(rmse0))]);end;
disp(['RMSE with INFER: '  num2str(mean(rmse))]);
disp(['RMSE with PoiR: '  num2str(mean(rmse2))]);
if dataSource==1;disp(['mean correct probability with TRUE: ' num2str(mean(probHit0(:)))]);end;
disp(['mean correct probability with INFER: ' num2str(mean(probHit(:)))]);
disp(['mean correct probability with PoiR: '  num2str(mean(probHit2(:)))]);

%% binary prediction: whether any event happens in an interval
intervalCand = [1, 2];
aucAll = [];
for i = 1:length(intervalCand)
    interval = intervalCand(i);
    figure;
    if dataSource==1
        model0.ratePar = ratePar;
        [inst0] = binPredict(model0, markedEventsAll, Tint, Tend, interval);
        [xx0,yy0, th0, auc(1)] = perfcurve(inst0(:,1),inst0(:,2),1);
        plot(xx0,yy0,'r');
    end;
    [inst1, feaTest] = binPredict(model, markedEventsAll, Tint, Tend, interval);
    [xx1,yy1, th1, auc(2)] = perfcurve(inst1(:,1),inst1(:,2),1);
    hold on;
    plot(xx1,yy1,'b');
    [ instTrain, feaTrain] = binPredict(model, markedEventsTrain, Tint, 0, interval);
    bb = glmfit(feaTrain, instTrain(:,1), 'binomial');
    poissMean = exp(poissWeight(1)+feaTest*poissWeight(2:end));
    pp = 1- exp(-poissMean);
    [xx2,yy2, th2, auc(3)] = perfcurve(inst1(:,1), pp,1);
    plot(xx2,yy2,'k');
    if dataSource==1; legend('truth','inferred','poissReg');else legend('inferred','poissReg');end;
    disp(['AUC: TRUE vs. INFER vs. PoiR ' num2str(auc)]);
    aucAll = [aucAll; auc];
end;
