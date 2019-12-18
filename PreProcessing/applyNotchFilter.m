function outSig = applyNotchFilter(sig, Fs, frequency, bandFrequency)
% Apply a notch filter with a FFT transfrom around the frequecny and its
% harmonics.
% sig: input signal with columns = time and rows = electrodes;
% Fs: sampling rate
% frequency: frequency to be removed,
% bandFrequency: delta around the frequency

if nargin<=1, frequency=50; end;
if  nargin<=2, bandFrequency = 1; end


removeFrequency = [frequency 2*frequency 3*frequency];  %#notch frequency
fn = Fs/2;              %#Nyquist frequency
freqRatio = removeFrequency/fn;      %#ratio of notch freq. to Nyquist freq.


[nChannel,nSamp] = size(sig);
f = (1:nSamp)*Fs/nSamp;


for iFreq = 1:size(removeFrequency,2)
    mask(f>removeFrequency(iFreq)-iFreq*bandFrequency & f<removeFrequency(iFreq)+iFreq*bandFrequency) = 1;
    mask(f>(Fs-removeFrequency(iFreq))-iFreq*bandFrequency & f<(Fs-removeFrequency(iFreq))+iFreq*bandFrequency) = 1;
end

outSig = zeros(nChannel,nSamp);
for iChannel=1:nChannel
    fftx = fft(sig(iChannel,:));
    fftx(mask==1)=0;
    outSig(iChannel,:) = real(ifft(fftx));  
end

end
