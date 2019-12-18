function outSig = applyNotchFilter(obj,frequency, bandFrequency)

% Apply a notch filter with a FFT transfrom around the frequecny and its
% harmonics.
% frequency: frequency to be removed, 
% deltaFrequency: delta around the frequency

if nargin<=1, frequency=50; end;
if  nargin<=2, bandFrequency = 1; end


Fs = obj.samp;
removeFrequency = [frequency 2*frequency 3*frequency]; 
nSamp= obj.nSamp;
f = (1:nSamp)*Fs/nSamp;

for iFreq = 1:size(removeFrequency,2)
mask(f>removeFrequency(iFreq)-iFreq*bandFrequency & f<removeFrequency(iFreq)+iFreq*bandFrequency) = 1;
mask(f>(Fs-removeFrequency(iFreq))-iFreq*bandFrequency & f<(Fs-removeFrequency(iFreq))+iFreq*bandFrequency) = 1;
end
outSig = zeros(obj.nChannel,obj.nSamp);
for iChannel=1:obj.nChannel
fftx = fft(obj.pProcessedSignal(iChannel,:));
fftx(mask==1)=0;
outSig(iChannel,:) = real(ifft(fftx));

end


%% Output
if nargout==0
    obj.pProcessedSignal = outSig;
    obj.filteringDetails{size(obj.filteringDetails,2)+1} =  [ 'Notch filter at ' num2str(frequency) ];
    % Plot
    if (obj.pdisplay==1)
        plot(obj);
    end
else
    if (obj.pdisplay==1)
        plot(outSig);
    end
end

end