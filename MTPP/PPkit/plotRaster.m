function [  ] = plotRaster(events, fig)
%PLOTRASTER Summary of this function goes here
figure (fig);
hold on;        
    yLimits = ylim;
                yAxisMax = yLimits(2);
        yAxisMin = -0.03 * yAxisMax;
        ylim([yAxisMin, yAxisMax]);
        line(xlim, [0 0], 'Color', 'black'); 
        xRug(events(:,1), .02, 'red');
end

