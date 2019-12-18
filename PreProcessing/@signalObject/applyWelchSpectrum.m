function [powerSpectrum,frequencyVector,localSpectrum] = applyWelchSpectrum(obj, varargin)

% Compute a spectrum based on Welch algorithm

% Input
nPoint = obj.nSamp;
nChannel = obj.nChannel;
samplingFrequency = obj.samp;
p=inputParser;
addParameter(p,'nPointFFT',4096);
addParameter(p,'preprocessing','no');
addParameter(p,'filter',[],@(x) (size(x,2)==2));
addParameter(p,'nOverlap',0.5);
parse(p,varargin{:});
nPointFFT = p.Results.nPointFFT;
preprocessing = p.Results.preprocessing;
nOverlap = p.Results.nOverlap;

if ~isempty(p.Results.filter), preprocessing='filter';end

% main loop
powerSpectrum = zeros(nChannel,nPointFFT);
frequencyVector = samplingFrequency*(0:nPointFFT-1)/nPointFFT;

for iChannel = 1:nChannel
   
    sig = obj.processedSignal(iChannel,:);
    
    %% Preprocessing
    sig = sig-mean(sig);
    
    if (strcmp(preprocessing,'filter') && ~isempty(filterCutoffFreq))     
        %disp(['Filter : [' num2str(filterCutoffFreq) '] Hz']);
        [b,a] = butter(4,filterCutoffFreq/ (samplingFrequency/2) );
        sig  = filter(b,a,sig);
    end
    
    if (strcmp(preprocessing,'everett') )  % for bipolar only, should be 0 for unipolar data
        % filter 40-250
        if obj.samp>500,
            [b,a] = butter(4,[40 250]/ (samplingFrequency/2) );
        else 
            [b,a] = butter(4,[40 ]/ (samplingFrequency/2) ,'high');
        end
        y  = filter(b,a,sig);
        % rectification
        sig=abs(y);
        
        % low pass 20Hz
        [b,a] = butter(2,20/ (samplingFrequency/2) );
        sig  = filter(b,a,sig);
    end
    
    
    %% Welch windowing and FFT computation
    
    powerSpectrumLocal = zeros(1,nPointFFT);
    nWindow = 0;
    iStartPoint = 1;
    endPoint = nPoint-nPointFFT*0.8;        % the last window contains at least 80% of signal, that will be padded with 0
    if (endPoint)<0, endPoint=1;end;
    iCp = 1;
        
    while iStartPoint<=endPoint
        spectrum = computeFFT(sig,samplingFrequency,nPointFFT,iStartPoint);
        if iChannel==1
            localSpectrum{iCp}.spectrum = zeros(nChannel,nPointFFT);
        end
        
        localSpectrum{iCp}.spectrum(iChannel,:) = spectrum;
        localSpectrum{iCp}.start = iStartPoint;
        localSpectrum{iCp}.stop = iStartPoint+nPointFFT;
        powerSpectrumLocal = powerSpectrumLocal + spectrum;
        nWindow = nWindow+1;
        iStartPoint = iStartPoint+nPointFFT*nOverlap;
        iCp=iCp+1;
    end    
    powerSpectrum(iChannel,:) = powerSpectrumLocal/nWindow;
    
end

function [powerSpectrum,frequencyVector,S] = computeFFT(sigIn,samplingFrequency,nPointFFT,startingPoint)


%% input and initialization
if nargin<=1, samplingFrequency = 1; end;
if nargin<=2, nPointFFT = size(sigIn,2); end;
if nargin==4, sigIn = sigIn(:,startingPoint:end); end;
hannWindow = 1;
nSig = size(sigIn,2);

nbSig = size(sigIn,1);
%disp('FFT ...');


for iSig = 1:nbSig
    sig = sigIn(iSig,:);
    sig=sig-mean(sig);
%% Zero Padding if needed
if (nSig<nPointFFT)
    sig2 = zeros(1,nPointFFT);
    sig2(1:nSig) = sig;
    sig = sig2;
elseif (nSig>nPointFFT)
    sig = sig(1:nPointFFT);
end


%% Hanning window if required
if hannWindow
    W = hanning(nPointFFT);
    sig = sig.*W';
end

%% Spectrum computation
S = fft(sig,nPointFFT); % Spectrum
frequencyVector = samplingFrequency*(0:nPointFFT-1)/nPointFFT;
powerSpectrum(iSig,:) = S.*conj(S)/nPointFFT;
end
