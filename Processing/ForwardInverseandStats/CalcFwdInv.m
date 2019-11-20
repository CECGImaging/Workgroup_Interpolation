function [fwd_bsp,inv_egm,fwdStats,invStats,fwdOutputs,invOutputs] = CalcFwdInv(torso,heart,EGM,BSP,expNotebook)
%Handles Calculating the Forward and Inverse depending on provided
%parameters
%Inputs:
% torso -Torso geometry with .node points and .face faces. Assumes trisurf
% heart -Heart geometry with .node points and .face faces
% EGM -Ground truth EGM data, provide if inverse statistics are to be calculated
% BSP -Ground truth BSP data, provide if forward statistics are to be calculated
% 
%expNotebook -struct array with parameters for the forward, inverse, and stats processes with fields:
%     notebook.verbose      -If the function should print things
%     notebook.fwdMethod    -Forward method to be used. 
%                               -1: skips forward function and looks for 
%                                provided forward matrix as 
%                                notebook.fwdSettings.fwd_mat
%                                 0: user provided forward function.
%                                 1: SCIRun BEM forward
%     notebook.invMethod    -Inverse Method to be used
%                                 -1: skips inverse calculation. Sets invEGM 
%                                     as all zeros
%                                 0: user provided inverse function
%                                 1: jaume_coll_font tikhonov implementation
%                                 https://github.com/CECGImaging/inverseecg
%     notebook.fwdFunct-if user defined
%  [fwd_bsp,fwdOutputs] = expNotebook.fwdFunct(torso,heart,EGM,fwdSettings);
% 
%     notebook.invFunct-if user defined
%     [inv_egm,invOutputs] = expNotebook.invFunct(fwd_mat,BSP,invSettings);
% 
%     notebook.doStats -true or false, detertmins if stats are to be
%                       calculated
%                             -Calculates RMSE, Temporal correlation, and spatial correlation
%     notebook.statFunct -if user defined, calculates extra stats
%             extraStats = expNotebook.statFunct(BSP,fwd_bsp,statSettings);
%             extraStats = expNotebook.statFunct(EGM,inv_egm,statSettings);
%     notebook.statSettings-passed to user defined stats
% 
%     fwdSettings.fwd_mat%if provided
%     fwdSettings %the rest of this struct depends on the fwd function used
% 
%     invSettings -depends on desired inv method.
%     Settuings to consider for default method 1:
%         invSettings.regMat%if specific regulatization desired
%         invSettings.weightMat%if weight mat desired, default identity
%         invSettings.lcurveParams
%         invSettings.underdetermined
%         invSettings.frobenius

if ~isfield(expNotebook,'verbose')||isempty(expNotebook.verbose), verb = true; else, verb = expNotebook.verbose; end
if verb
    fprintf('FWD INV Manager startup. Checking inputs\n');
end

if ~isfield(expNotebook,'doStats')||isempty(expNotebook.doStats), expNotebook.doStats = true; end

invSettings = expNotebook.invSettings;
fwdSettings = expNotebook.fwdSettings;

if isfield(torso,'pts')
    torso.node = torso.pts;
    torso = rmfield(torso,'pts');
end
if isfield(heart,'pts')
    heart.node = heart.pts;
    heart = rmfield(heart,'pts');
end
if isfield(torso,'fac')
    torso.face = torso.fac;
    torso = rmfield(torso,'fac');
end
if isfield(heart,'fac')
    heart.face = heart.fac;
    heart = rmfield(heart,'fac');
end



%calc fwd
forward = true;

switch expNotebook.fwdMethod
    case -1
        
        if verb
            fprintf('Skipping Forward Calculation\n');
        end

        try
            fwdOutputs.fwd_mat = fwdSettings.fwd_mat;
            fwd_bsp = fwdOutputs.fwd_mat*EGM;
            fprintf('Using Provided FWD Matrix\n');
        catch
            fwd_bsp = zeros(size(BSP));
            forward = false;
            fprintf('FWD solution disabled but no FWD matrix provided in fwdSettings.fwd_mat.\nSubsequent steps likely to fail.');
        end
        
    case 0%User proivided method
        if verb
            fprintf('User Defined Forward Calculation\n');
        end
        [fwd_bsp,fwdOutputs] = expNotebook.fwdFunct(torso,heart,EGM,fwdSettings);
    case 1%Standard BEM using scirum matlab BEM
        if verb
            fprintf('BEM Forward Calculation\n');
        end
        [fwd_bsp,fwd_mat,opt] = Handle_BEM_FWD(torso,heart,EGM,fwdSettings);
        fwdOutputs.fwd_mat = fwd_mat;
        fnames = fieldnames(opt);
        for f = 1:length(fnames)
            fwdOutputs.(fnames{f}) = opt.(fnames{f});
        end
end


%processFwd
if isfield(fwdSettings,'fwdPostProcess')
    if verb
        fprintf('Post Processing FWD Results using user defined function\n');
    end
    [fwd_bsp,fwdOutputs.fwd_mat] = fwdSettings.fwdPostProcess(fwd_bsp,BSP,fwdOutputs.fwd_mat,fwdSettings);
    
end


%calc inverse
inverse = true;
switch expNotebook.invMethod
    case -1%Do not do an inverse
        if verb
            fprintf('Skipping Inverse Calculation\n');
        end
        inverse = false;
        inv_egm = zeros(size(EGM));
        invOutputs.note = 'skipped inv calculation';
    case 0%user provided method
        if verb
            fprintf('User Defined Inverse Calculation\n');
        end
        [inv_egm,invOutputs] = expNotebook.invFunct(fwd_mat,BSP,invSettings);
    case 1 %Jaume tikhonov
        
        if verb
            fprintf('Jaume Inverse Calculation\n');
        end
        if isfield(invSettings,'regMat')
            regMat = invSettings.regMat;
        else
            if verb
                fprintf('No Reg matrix given, using laplacian.\n');
            end
            heartAdjacency = computeAdjacencyMatrix(heart,2);
            weightFunct = @(index) heartAdjacency(index,:);
            [D,H] = meshVolDiffHessMatrix(heart,weightFunct);
            heartLaplacian = LaplacianMatrixFromHessianMatrix(H);
            regMat = heartLaplacian;
        end
        if ~isfield(invSettings,'weightMat')||isempty(invSettings.weightMat)
            invSettings.weightMat = 0;
        end
        weightMat = invSettings.weightMat;
        if weightMat == 0%If 0 selected then generate identity weight matrix
            weightMat = eye(size(fwdOutputs.fwd_mat,1));
        end
        lambdaVals = 10.^(2*linspace(invSettings.lcurveParams(1),invSettings.lcurveParams(2),invSettings.lcurveParams(3)));
        determined = invSettings.underdetermined;
        if ~isfield(invSettings,'frobenius')
            frobenius = true;
        else
            frobenius = invSettings.frobenius;
        end
        [inv_egm, lambdaCornerIX, ~, rho, eta, kappa] = tikhonov_jcf(fwdOutputs.fwd_mat,regMat,weightMat,BSP,lambdaVals,determined,frobenius);
        invOutputs.lambdaCornerIX = lambdaCornerIX;
        invOutputs.rho = rho;
        invOutputs.eta = eta;
        invOutputs.kappa = kappa;
end



%calc stats/metrics
if expNotebook.doStats
    if verb
        fprintf('Calculating Statistics from Results');
    end
    %fwd stats
    if forward
        if verb
            fprintf(' . ');
        end
        %signed difference
        fwd_diff = BSP-fwd_bsp;
        %RMSE
        RMSE = sqrt(sum(sum(fwd_diff.^2))/length(BSP(:)));
        %PCC vals for each time point
        PCCs = zeros(1,size(BSP,2));
        for t = 1:size(BSP,2)
            PCCs(:,t) = corr2(BSP(:,t),fwd_bsp(:,t));
        end
        %mean PCC val for all time points
        PCC_mean = mean(PCCs);
        
        %cross correlation (normalized) for each electrode
        corrVals = zeros(size(BSP,1),1);
        for e = 1:size(BSP,1)
            corrVals(e,:) = xcorr(BSP(e,:),fwd_bsp(e,:),0,'coeff'); 
        end
        %mean cross correlation across electrodes
        corr_mean = mean(corrVals);
        %extra used defined stats functions
        if isfield(expNotebook,'statFunct')
            statSettings = expNotebook.statSettings;
            extraStats = expNotebook.statFunct(BSP,fwd_bsp,statSettings);
        else
            extraStats = [];
        end
        if verb
            fprintf(' . ');
        end
        

        


    else
        fwd_diff = [];
        RMSE = [];
        PCCs = [];
        PCC_mean = [];
        corrVals = [];
        corr_mean = [];
        extraStats = [];
    end
    if verb
        fprintf(' . ');
    end
    fwdStats.diff = fwd_diff;
    fwdStats.RMSE = RMSE;
    fwdStats.PCCs = PCCs;
    fwdStats.PCC_mean = PCC_mean;
    fwdStats.corrVals = corrVals;
    fwdStats.corr_mean = corr_mean;
    fwdStats.extraStats = extraStats;

    %inverse
    if inverse
        if verb
            fprintf(' . ');
        end
        %signed difference
        inv_diff = EGM-inv_egm;
        %RMSE
        RMSE = sqrt(sum(sum(inv_diff.^2))/length(EGM(:)));
        %PCC vals for each time point
        PCCs = zeros(1,size(EGM,2));
        for t = 1:size(EGM,2)
            PCCs(:,t) = corr2(EGM(:,t),inv_egm(:,t));
        end
        %mean PCC val for all time points
        PCC_mean = mean(PCCs);
        
        %cross correlation (normalized) for each electrode
        corrVals = zeros(size(EGM,1),1);
        for e = 1:size(EGM,1)
            corrVals(e,:) = xcorr(EGM(e,:),inv_egm(e,:),0,'coeff'); 
        end
        %mean cross correlation across electrodes
        corr_mean = mean(corrVals);
        %extra used defined stats functions
        if isfield(expNotebook,'statFunct')
            statSettings = expNotebook.statSettings;
            extraStats = expNotebook.statFunct(EGM,inv_egm,statSettings);
        else
            extraStats = [];
        end
        if verb
            fprintf(' . ');
        end



    else
        inv_diff = [];
        RMSE = [];
        PCCs = [];
        PCC_mean = [];
        corrVals = [];
        corr_mean = [];
        extraStats = [];
    end
    invStats.diff = inv_diff;
    invStats.RMSE = RMSE;
    invStats.PCCs = PCCs;
    invStats.PCC_mean = PCC_mean;
    invStats.corrVals = corrVals;
    invStats.corr_mean = corr_mean;
    invStats.extraStats = extraStats;
else
    fwdStats = [];
    invStats = [];
    if verb
        fprintf('Skipping stats calculation');
    end

end
if verb
    fprintf(' .\nDone\n');
end

end



