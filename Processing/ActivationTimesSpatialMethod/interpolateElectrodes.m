function [ newPots ] = interpolateElectrodes( geomObj, pots )
%INTERPOLATEELECTRODES Replace all NaN electrodes by interpolations
[nmbElec,nmbTime] = size(pots);
newPots = pots;

% Create the neighbor matrix
neighborMatrix = zeros(geomObj.NrVertices);
for i=1:geomObj.NrVertices;
        ind = [i; getNeighbors(i,geomObj.faces)];
        neighborMatrix(i,ind)=1;
        neighborMatrix(ind,i)=1;
end

parfor i=1:nmbElec
    for j=1:nmbTime
        if(isnan(pots(i,j)))
            % interpolate
            indNeighbors=findNonNANneighbors(pots(:,j),neighborMatrix,i,i);
            
%             trisurf(geomObj.faces,geomObj.vertices(:,1),geomObj.vertices(:,2),geomObj.vertices(:,3),'FaceColor',[.49 1 .63],'Marker','o','MarkerSize',1,'MarkerEdgeColor','k','MarkerFaceColor','r','FaceAlpha',0.8); axis equal; title('New surface');
%                 hold on; scatter3(geomObj.vertices(indNeighbors,1),geomObj.vertices(indNeighbors,2),geomObj.vertices(indNeighbors,3)); hold on;
%     scatter3(geomObj.vertices(i,1),geomObj.vertices(i,2),geomObj.vertices(i,3),'o','MarkerFaceColor','r');
%             
            newPots(i,j) = interpolateElec(geomObj.vertices(i,:), geomObj.vertices(indNeighbors,:), pots(indNeighbors,j));
        end
    end
end

end

function [indNeighbors ignoreVertices] = findNonNANneighbors(potsAtTime,neighborMatrix,currentVertex,ignoreVertices)
    indNeighborsTemp=find(neighborMatrix(currentVertex,:)==1);
    indNeighborsTemp = setdiff(indNeighborsTemp,ignoreVertices);
    indNeighbors = [];
    % First, increase ignoreVertices list to all neighbours that are NaN
    for k=1:length(indNeighborsTemp)
        if(isnan(potsAtTime(indNeighborsTemp(k))))
            ignoreVertices = unique([ignoreVertices indNeighborsTemp(k)]);
        end
    end
    % Then, look for their neighbours until a non-NaN is found
    for k=1:length(indNeighborsTemp)
        if(isnan(potsAtTime(indNeighborsTemp(k))))
            [indNeighborsNew ignoreVertices] = findNonNANneighbors(potsAtTime,neighborMatrix,indNeighborsTemp(k),ignoreVertices);
            indNeighbors = unique([indNeighbors indNeighborsNew]);
        else
            indNeighbors = [indNeighbors indNeighborsTemp(k)];
        end
    end
end

function value = interpolateElec(currentPoint, neighborPoints, neighborPots)
    for i = 1:size(neighborPoints,1)
        neighborDistances(i) = norm(neighborPoints(i,:) - currentPoint);
    end
    % Addition 2019: only keep neighbors until a max distance is found
    if true % set to true if you want to have a max distance
        maxDistance = 15; % in millimeter (or whathever your units are)
        closenodes = find(neighborDistances<=maxDistance);
%         disp(['Keeping ' num2str(length(closenodes)) ' out of ' num2str(length(neighborDistances)) ' nodes for interpolation']);
        neighborPoints = neighborPoints(closenodes,:);
        neighborPots = neighborPots(closenodes,:);
        neighborDistances = neighborDistances(closenodes);
    end
    totalDistance = sum(neighborDistances);
    weights = (1-(neighborDistances/totalDistance));
    weights = weights./sum(weights);
    if length(neighborDistances>0)
        value = sum(neighborPots.*weights');
    else
        value = NaN;
    end
%     disp([num2str(value) ' with ' num2str(length(neighborDistances)) 'neighbors'] );
%     value = sum(neighborPots)/length(neighborPots);
end