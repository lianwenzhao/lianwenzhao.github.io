function [ fig ] = plotPiecewiseRate(rateTrue, fig, color)
%PLOTPIECEWISERATE Summary of this function goes here
minDelta = 0;
if nargin<3
    color ='r';
end;
if nargin<2
    fig = figure;
end;
figure(fig);
hold on;
i=2;
if size(rateTrue, 2) >2
    while (i<size(rateTrue,1))
        if (rateTrue(i,1)-rateTrue(i-1,1))>minDelta
            shadedErrorBar([rateTrue(i-1,1), rateTrue(i,1)], [rateTrue(i,2), rateTrue(i,2)],  [rateTrue(i,3), rateTrue(i,3)], color);
            plot([rateTrue(i,1), rateTrue(i,1)], [rateTrue(i,2), rateTrue(i+1,2)],[color '--'],'lineWidth',1);
        end;
        i = i+1;
    end;
    shadedErrorBar([rateTrue(i-1,1), rateTrue(i,1)], [rateTrue(i,2), rateTrue(i,2)],  [rateTrue(i,3), rateTrue(i,3)], color);
    xylim = axis();
    axis([min(xylim(1), min(rateTrue(:,1))), max(max(rateTrue(:,1))), 0, max([max(rateTrue(:,2))*1.1, 1.5, xylim(4)])]);
    drawnow;
else
    while (i<size(rateTrue,1))
        if (rateTrue(i,1)-rateTrue(i-1,1))>minDelta
            plot([rateTrue(i-1,1), rateTrue(i,1)], [rateTrue(i,2), rateTrue(i,2)],color,'lineWidth',1);
            plot([rateTrue(i,1), rateTrue(i,1)], [rateTrue(i,2), rateTrue(i+1,2)],[color '--'],'lineWidth',1);
        end;
        i = i+1;
    end;
    plot([rateTrue(i-1,1), rateTrue(i,1)], [rateTrue(i,2), rateTrue(i,2)],color,'lineWidth',1);
    xylim = axis();
    axis([min(xylim(1), min(rateTrue(:,1))), max(max(rateTrue(:,1))), 0, max([max(rateTrue(:,2))*1.1, 1.5, xylim(4)])]);
    drawnow;
end

