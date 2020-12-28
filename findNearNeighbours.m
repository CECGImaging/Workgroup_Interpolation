function [neighbour, distance] = findNearNeighbours(tri,whichPoint)

distance = [];
neighbour = [];
% Manage multiple request
if isempty(whichPoint)
    neighbour=[];
    distance = 0;
    return;
end

if size(whichPoint,2)>1
    neighbour = [];
    for i=1:size(whichPoint,2)
        neighbour= [neighbour FindNeighbour(tri,whichPoint(i))];    %#ok
        neighbour = sort(unique(neighbour));
    end
    distance=[];
else
    face = vertexAttachments(tri,whichPoint);
    vertex = tri.Triangulation(face{1},:);
    if size(vertex,1)==1
        vertex = [vertex;vertex];
    end
    neighbour = intersect(sort(vertex(vertex ~= whichPoint)), [1:size(tri.X,1)]);
    for iN =1:size( neighbour,1)
        distance(iN) = norm(tri.X(whichPoint,:)-tri.X(neighbour(iN),:));
    end
    
end
if isempty(distance) && isempty(neighbour)
    distance = 0;
    neighbour = 0;
end

end