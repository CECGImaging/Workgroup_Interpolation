function [ timeInMs, timeInMs_Unsmoothed, indexWithEarliestTime, indexWithEarliestTime_Unsmoothed, timeconfidence, maxSlope ] = getTimeFromElectrogram( pots, timewindow, samplingrate, settings, postsmoothingsettings  )
%GETTIMEFROMELECTROGRAM Summary of this function goes here
%   Estimations local depolarization (or, by giving potentials with a negative sign,
%   local repolarization) from a set of electrograms (and their possible
%   spatial relation)
% Input variables:
%   pots: the set of electrograms (each row is one electrogram)
%   timewindow: the time window on which to compute the activation time
%   samplingrate: what it says
%   settings: method-specific settings, such as
%   settings.method: 'Spatiotemporal' ('Erem') or 'Temporal only'
%       settings.method = 'spatiotemp' || 'temponly';
%   for 'spatiotemp':
%       settings.Dtan = Dtan from meshsurfdiffhessmatrix()
%       settings.twindow = timewindow for spatiotemporal method (41 works
%       for 2048Hz, but see detailed studies on this)
%   for 'temponly':
%       settings.splineSmoothingParam = setting for smoothing for splines
%       per electrogram (0.001 works)
%   postsmoothingsettings: method-specific settings for post-estimation
%   spatial smoothing
%       postsmoothingsettings.dopostsmoothing = smoothing or no
%       smoothing (post-estimation)
%       postsmoothingsettings.Ltan = surface laplacian (e.g. from
%       mesh_laplacian())
%       postsmoothingsettings.regparamrange = range or single value for
%       smoothing (range seems not work nicely...)
%       postsmoothingsettings.geomobj = geom.Heart; necessary for
%       interpolating
%
% Output variables:
%   timeInMs: row-wise activation times (one time per electrogram). If
%   post-estimation spatial smoothing is used, these are the smoothed
%   times.
%   timeInMs_Unsmoothed: the unsmoothed activation times. (Identical to
%   TimeInMs if no post-estimation smoothing was applied.)
%   indexWithEarliestTime: row index of electrogram with earliest
%   time; based on smoothed times if post-estimation smoothing was used.
%   indexWithEarliestTime_Unsmoothed: location with earliest
%   time when no smoothing is applied (equal to previous variable if no
%   smoothing was used at all)
%   timeconfidence: confidence measure (Josselin method)

if ~exist('postsmoothingsettings','var'); postsmoothingsettings.dopostsmoothing = false; end;
if isempty(timewindow); timewindow = 1:size(pots,2); end;
nmbLeads = size(pots,1);
showSmoothingFigure = false;
maxSlope = NaN(nmbLeads,1);

% ESTIMATION OF ACTIVATION TIME
switch settings.method
    case 'spatiotemp'
        if showSmoothingFigure; figure; end;
        warning('Doing some smoothing for spatiotemp approach! Smoothing param = used defined');
        % some smoothing for spatiotemp
        smoothedPots = NaN(size(pots));
        parfor i = 1:size(pots,1)
            if ~any(isnan(pots(i,timewindow)));
                f = fit(timewindow',pots(i,timewindow)','smoothingspline','SmoothingParam', settings.splineSmoothingParam);
                smoothedPots(i,timewindow) = feval(f,timewindow)';
                smoothedDiffed = differentiate(f,timewindow); 
                maxSlope(i) = max(abs(smoothedDiffed));
                if showSmoothingFigure; plot(pots(i,timewindow),'r'); hold on; plot(smoothedPots(i,timewindow),'b'); plot(smoothedDiffed(i,:),':b'); title(['Max slope = ' num2str(slope)]); hold off; pause; end;
            end
        end
        [acttimesunsmoothedIndex,~,~,~,timeconfidence] = spatiotemporalActtimes(smoothedPots(:,timewindow),settings.Dtan,settings.twindow);
        % original (nonsmoothed) version: [acttimesunsmoothedIndex,~,~,~,timeconfidence] = spatiotemporalActtimes(pots(:,timewindow),settings.Dtan,settings.twindow);
        acttimesunsmoothedIndex = acttimesunsmoothedIndex+timewindow(1);
    case 'temponly'
        acttimesunsmoothedIndex = NaN(nmbLeads,1);
        timeconfidence = acttimesunsmoothedIndex;
        smoothedPots = NaN(size(pots));
        smoothedDiffedMin = NaN(nmbLeads,1);
%         thresholdValue = 0.250* nanmean(nanstd(pots(:,timewindow),0,2)); % more or less trial and error setting...
        for i = 1:1:size(acttimesunsmoothedIndex,1)
            if ~any(isnan(pots(i,timewindow)));
%                 thresholdValue = 0.0050* nanstd(pots(i,:)); % threshold define over the WHOLE signal (QRS + T wave) of THIS node
                thresholdValue = 1.0* nanstd(pots(i,:)); % threshold define over the WHOLE signal (QRS + T wave) of THIS node
                f = fit(timewindow',pots(i,timewindow)','smoothingspline','SmoothingParam', settings.splineSmoothingParam);
                smoothedDiffed = differentiate(f,timewindow);
                maxSlope(i) = max(abs(smoothedDiffed));
                smoothedPots(i,timewindow) = feval(f,timewindow)';
                if max(abs(smoothedPots(i,timewindow))) > thresholdValue % only take leads into account with a substantial potential (> threshold)
                    [blah,depind] = min(smoothedDiffed);
                    acttimesunsmoothedIndex(i) = timewindow(1)+depind;
                    belowzeroindices = find(smoothedDiffed<0);
                    auc = sum(smoothedDiffed(belowzeroindices));
                    timeconfidence(i) = blah/auc;
                end
            end
        end
end

% From index to time-in-ms:
timeInMs_Unsmoothed = acttimesunsmoothedIndex/samplingrate * 1000;
[~ , indexWithEarliestTime_Unsmoothed] = min(timeInMs_Unsmoothed);

% IF REQUIRED, POST-ESTIMATION SPATIAL SMOOTHING
if postsmoothingsettings.dopostsmoothing
    % Check whether there are NaNs and fill them:
    timeInMs = timeInMs_Unsmoothed;
    if any(isnan(timeInMs)); timeInMs = interpolateElectrodes(postsmoothingsettings.geomobj, timeInMs); end;
    timeInMs=smoothactivationtimes(postsmoothingsettings.Ltan,timeInMs,postsmoothingsettings.regparamrange);
    [~ , indexWithEarliestTime] = min(timeInMs);
else
    % if no post-estimation smooting: unsmoothed and 'smoothed' times are
    % the same
    timeInMs = timeInMs_Unsmoothed;
    indexWithEarliestTime = indexWithEarliestTime_Unsmoothed;
end

end

