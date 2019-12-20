function outSigOut=applyWaveletFilter(obj, nWavelet, filterband)

% Apply wavelet bandpass filtering
% the decompostion is performed over 'nWavelet' (default=10); and the
% filtering within 'filteringBand' (default=[3 12])

% input
if (nargin<=1 | isempty(nWavelet)), nWavelet=10;end
if nargin<=2,filterband=[3 12];end;

% local variable for multicore processing
outSigOut = zeros(obj.nChannel,obj.nSamp);
sigTemp = obj.pProcessedSignal;
nSamp  = obj.nSamp;
samp=obj.samp;

% main loop
parfor iChannel = 1:obj.nChannel
    sig = sigTemp(iChannel,:);
    [C,L]= wavedec(sig,nWavelet,'coif4');
    outSig=zeros(1,nSamp);
    for iWavelet=1:nWavelet,
        d= wrcoef('d',C,L,'coif4',iWavelet);
        dfft=d.*hanning(nSamp)';
        S = fft(dfft);
        f = samp*(0:nSamp-1)/nSamp;
        D = S.*conj(S)/nSamp;
        if max(D(f<=filterband(2) & f>=filterband(1)))== max(D)     %AF 3-12 3-8
            outSig = outSig+d;
        end
    end;
    outSigOut(iChannel,:) = outSig;
end


% Output
if nargout==0
    obj.pProcessedSignal = outSigOut;
    obj.filteringDetails{size(obj.filteringDetails,2)+1} =  [ 'Wavelet filter coif4: N=' num2str(nWavelet) ', band: ' num2str(filterband) ];
    % Plot
    if (obj.pdisplay==1)
        plot(obj);
    end
else  
    if (obj.pdisplay==1)
        plot(outSigOut);
    end
end

end