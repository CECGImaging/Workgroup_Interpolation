function outSig = applySgolayFilter(obj,timeWindow,order)

% Apply a Savitzy-Golay filter of n-order over timeWindow ms


if nargin<=1, timeWindow=41; end;
if  nargin<=2, order = 3; end

outSig = sgolayfilt(obj.pProcessedSignal',order,timeWindow);
outSig = outSig';


%% Output
if nargout==0
    obj.pProcessedSignal = outSig;
    obj.filteringDetails{size(obj.filteringDetails,2)+1} =  [ 'Savitzy-Golay filter over ' num2str(timeWindow)/obj.samp*1000 ' ms' ];
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