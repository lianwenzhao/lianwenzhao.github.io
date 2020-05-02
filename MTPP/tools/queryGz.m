function [g, dg ] = queryGz( z, gz, dgz )
%QUERYGZ Summary of this function goes here
% assert( z>=-500, ['out of bound! z = ' num2str(z)]);
% g= gz(floor(-z/0.001)+1);
% dg = dgz(floor(-z/0.001)+1);

if z>-700
    g= gz(floor(-z/0.001)+1);
    dg = dgz(floor(-z/0.001)+1);
else
    g = gz(700/0.001+1) + dgz(700/0.001) * (z+700);
    dg = dgz(700/0.001);
end;