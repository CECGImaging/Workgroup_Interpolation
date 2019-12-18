%This is an example script to demonstrate the filters Laura has implemented and how to use
%the signalObject class

clc
close all
clear all


%load some data
load('Run0038-ts_1_1_1')
load('badTankNodeMap')

%% high-frequency removal
%Notch Filter Removal
Signal = signalObject(Sig,'samp',1000,'label','Run0038','badLeads',tankBadNodeMap);
Signal.applyNotchFilter(60,1); % Use 50 for Bordeaux Data
plot(Signal,5) % Can visualize difference between filtered and unfiltered signals

% Savitzy-Golay high-frequency removal
Signal.resetProcessedSignal; % Remove previous filter
Signal.applySgolayFilter(17,3); % Use 41,3 for Bordeaux Data
plot(Signal,5) % Can visualize difference between filtered and unfiltered signals


%% base-line drift removal
% Wavelet
Signal.resetProcessedSignal; % Remove previous filter
Signal.removeBaseline('wavelet'); % Uses default parameters
plot(Signal,5) % Can visualize difference between filtered and unfiltered signals

% Savitzy-Golay BDR
Signal.resetProcessedSignal;
Signal.removeBaseline('savitzkygolay'); % Uses default parameters
plot(Signal,5) % Can visualize difference between filtered and unfiltered signals
       
% Spline - could be improved by using pre-defined beats for isoelectric
% point. 
Signal.resetProcessedSignal;
Signal.removeBaseline('sinusspline',[100 320]); % 100 = ms prior to peaks to select iso-electric points, 320 is approx ms between each beat
plot(Signal,5) % Can visualize difference between filtered and unfiltered signals

% Simple BDR - I haven't implemented this into the signalObject yet. Need
% to load a beat segmentation file 
newSig = Sig;
for beat = 1:nbBeats
    if(beat<nbBeats)
        Signal = signalObject(Sig(:,beats{beat}.fromIdx:beats{beat+1}.fromIdx),'samp',1000,'badleads',badleads);
        Signal.minusConstant(mean(Signal.processedSignal(:,1:17)')); % 41 for Bordeaux
        newSig(:,beats{beat}.fromIdx:beats{beat+1}.fromIdx)=Signal.processedSignal;
    else
        Signal = signalObject(Sig(:,beats{beat}.fromIdx:end),'samp',samp,'label','pacing','badleads',badleads);
        Signal.minusConstant(mean(Signal.processedSignal(:,1:Notch)'));
        newSig(:,beats{beat}.fromIdx:end)=Signal.processedSignal;
    end
end
figure()
plot(Sig(5,:),'b');
hold on
plot(newSig(5,:),'k')

    
