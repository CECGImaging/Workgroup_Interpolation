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


grndTruth = load([dataDir,'EGMTrutBeats.mat']);



fNameBase = 'X_Tikh_Case2_Pacing_Int_5_';