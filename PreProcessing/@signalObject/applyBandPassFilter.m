function outsig = applyBandPassFilter(obj,filterCutoffFreq,N)
% Apply a FIR band pass filter

if nargin==2, N=200;end
%calculate filter coefficients
[b,a] = fir1(N,filterCutoffFreq/(obj.samp/2));
%apply filter
outsig= zeros(obj.nChannel,obj.nSamp);

parfor iChannel=1:obj.nChannel
outsig(iChannel,:)= filtfilt(b,a,obj.pProcessedSignal(iChannel,:));
end
if nargout==0
    obj.pProcessedSignal = outsig;
    obj.filteringDetails{size(obj.filteringDetails,2)+1} =  [ 'Band-pass filter [' num2str(filterCutoffFreq(1)) '-' num2str(filterCutoffFreq(2)) ']' ];
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