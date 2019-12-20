function [D,H,Ltan, Dstruct, Hstruct, Ltanstruct,AdjMtrxStruct] = getSpatialOperators(geomsub,maxdepth)
heart.node = geomsub.vertices';
heart.face = geomsub.faces';
[AdjMtrx,AdjMtrxStruct] = computeAdjacencyMatrixRecur(heart, maxdepth);
% compute gradient estimator
for depth=1:maxdepth
    disp(['Calculating operators for depth = ' num2str(depth) '/' num2str(maxdepth)]);
    AdjMtrxTemp =AdjMtrxStruct{depth};
    wghFcn = @(indx) AdjMtrxTemp(indx,:);
    [Dstruct{depth}, Hstruct{depth}] = meshVolDiffHessMatrix(heart,wghFcn);	
    Ltanstruct{depth}=LaplacianMatrixFromHessianMatrix(Hstruct{depth});
end
D = Dstruct{end}; H = Hstruct{end}; Ltan = Ltanstruct{end};