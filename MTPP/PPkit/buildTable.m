function [ piecePop, pieceAll, piecePopRep ] = buildTable( pieceAll )
%BUILDTABLE Summary of this function goes here
%%%input: feature set from all individuals
%%%output:
%%%piecePop: feature set for groups (population level)
%%% pieceHt: hash table, key is the string for feature, value is latent
%%% function value
tmp = [];
for u=1:length(pieceAll)
    pieceAll{u}.hashKey = cell(size(pieceAll{u}.feature,1),1);
    tmp = [tmp; pieceAll{u}.feature];
    for i=1:length(pieceAll{u}.hashKey)
        pieceAll{u}.hashKey{i} = num2str(pieceAll{u}.feature(i,:));
    end;
end
piecePopRep = tmp;
piecePop = unique(tmp, 'rows');


% a=num2str([41;4])
% b=num2str([4;14])
% c = containers.Map;
% c(a) = 1
% c(b) = 2
% keys(c)
% values(c)

