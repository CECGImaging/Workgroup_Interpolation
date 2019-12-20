%%
clear
close all
clc
%%
%Script to compare the inverse solutions to the the ghround truth and save
%out statistics for visualization


%Data directory where inverse results are stored
dataDir = '/home/sci/jbergquist/projects/WorkingData/InterpolationGroup/TikhonovRegResults/';

%Output directory to save visualization data to
outDir = '/home/sci/jbergquist/projects/WorkingData/InterpolationGroup/TikhonovRegResults/Vis/';


grndTruth = load([dataDir,'EGMTrueBeats.mat']);



fNameBase = 'X_Tikh_Case2_Pacing_Int_5_';

for pt2 = 1:5
    for interp = 1:7
        
        try
           load([dataDir,fNameBase,sprintf('%d_%d.mat',pt2,interp)]);
           
           for b = 1:31
               currBeat = X_Tikh{b};
               
               [ RMSE, RMSE_space, RMSE_time, PCC, mean_PCC, CorrVal, mean_CorrVal,Diff,maxDiff,minDiff ] = CalculateSignalStatistics( currBeat.VeSingLam(1:108,:), grndTruth.EGMtrue{b} );
               save([outDir,fNameBase,sprintf('%d_%d_RMSE_Space_b%d.mat',pt2,interp,b)],'RMSE_space');
               save([outDir,fNameBase,sprintf('%d_%d_RMSE_CorrVal_b%d.mat',pt2,interp,b)],'RMSE_space');
           end
            
        catch
            fprintf(['Failed at ',fNameBase,sprintf('%d_%d\n',pt2,interp)]);
        end
        
        
    end
end