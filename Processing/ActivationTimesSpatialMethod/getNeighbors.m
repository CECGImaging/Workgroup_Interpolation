function [ neighbors ] = getNeighbors( vertex, TRI )
%GETNEIGHBORS Returns the neighbors of vertex in the triangularizaiton TRI
[allEdges,jj]=sortrows(sort([TRI(:,[1 2]); TRI(:,[1 3]); TRI(:,[2 3])],2));
connectedPoints = unique(allEdges,'rows');
% if point i is connected with point j then j is connected with i as well
connectedPoints = [connectedPoints; fliplr(connectedPoints)];
neighbors = connectedPoints(find(ismember(connectedPoints(:,1),vertex)),2);

end

