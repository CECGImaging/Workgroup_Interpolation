%% HELP:
%
%	This functionl computes the adjacency matrix on the graph defined in
%	the geometry "geom".
%	The adjacency matrix can be extended to all the neighbors of a node
%	that are at a distance "pathLength" from that node.
%
%
%	INPUT:
%		- geom - struct - input geometry. Assumed to be a triangular
%		surface mesh. The struc is composed of two fields: node and face.
%		- pathLength - int - number of jumps on the graph that define the
%		adjacency.
%
%	OUTPUT:
%		- AdjMtrx - <M,M>bool - matrix defining which nodes are neighbors
%		of which.
%
%	AUTHOR:
%		Jaume Coll-Font <jcollfont@gmail.com>
%   ADDAPTATION by Matthijs Cluitmans 9 dec 2016: returning all adjmatrices
%   for shorter paths as well

function [AdjMtrx,AdjMtrxStruct] = computeAdjacencyMatrixRecur(geom, pathLength)

	%% define
		[d,M] = size(geom.node);
		[~,numFac] = size(geom.face);
			
		if d ~= 3
			numFac = d;
			geom.face = geom.face';
			geom.node = geom.node';
		end
			
	%% create 1rst order adjacency matrix
		adjMtrx = false(M);
		for ii = 1:numFac;
			adjMtrx(geom.face(1,ii),geom.face(2,ii)) = 1;%distance(geom.face(1,ii),geom.face(2,ii));
			adjMtrx(geom.face(2,ii),geom.face(3,ii)) = 1;%distance(geom.face(2,ii),geom.face(3,ii));
			adjMtrx(geom.face(3,ii),geom.face(1,ii)) = 1;%distance(geom.face(3,ii),geom.face(1,ii));
			adjMtrx(geom.face(2,ii),geom.face(1,ii)) = 1;%distance(geom.face(2,ii),geom.face(1,ii));
			adjMtrx(geom.face(3,ii),geom.face(2,ii)) = 1;%distance(geom.face(3,ii),geom.face(2,ii));
			adjMtrx(geom.face(1,ii),geom.face(3,ii)) = 1;%distance(geom.face(1,ii),geom.face(3,ii));
        end

		
	%% extend adjacency matrix according to path length
	AdjMtrx = adjMtrx + eye(M);
    AdjMtrxStruct{1} = AdjMtrx;
	for rep = 1:(pathLength-1)
        disp(['Calculating adjacency matrix for depth = ' num2str(rep+1) '/' num2str(pathLength)]);
		for n =1:M
			neighbors = find(AdjMtrx(n,:));
			for m = neighbors
				AdjMtrx(n,:) = AdjMtrx(n,:)|adjMtrx(m,:);
			end
        end
        AdjMtrxStruct{rep+1} = AdjMtrx;
	end
end