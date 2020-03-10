function dist = EuclideanDistance(point1,point2)
% Calculates the Euclidean distance between two points
% Inputs point1 and point2 are (x,y,z) of two points

dist = sqrt( (point1(1)-point2(:,1)).^2 + (point1(2)-point2(:,2)).^2 + (point1(3)-point2(:,3)).^2 );
 
return
