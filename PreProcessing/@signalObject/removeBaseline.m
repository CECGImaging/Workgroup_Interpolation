function outSigOut  = removeBaseline(obj,method,parameter)

% Removes the baseline
%
% 4 methods can be used
% poly: use a polynomial interpolation, order of the polynome is
% proportionnal to the lenght of the signal: =2+1 order every 8 sec of
% signal
% input:
% signalObjet
% method: can be
% -filter: apply a bandpass filter between .5 150]Hz or sepcified as
% a third input
% -wavelet: apply wavelet filter between [.5 150]]Hz or sepcified as
% a third input
% -spline (use for VF): look for anchor points every 250ms (or sepcified as
% a third input in ms) and compute a spline interpolation
%
% if no output are asked, results are saved in the object


if (nargin==1), method = 'poly';end
if (nargin==2), parameter = []; end
outSigOut = zeros(obj.nChannel,obj.nSamp);

switch lower(method)
    case 'savitzkygolay'
        if isempty(parameter); parameter=[3 round(3000/(obj.samp/1000))]; end
        if(mod(parameter(2),2)==0);parameter(2) = parameter(2)+1;end 
        
        baseline = sgolayfilt(obj.processedSignal',parameter(1),parameter(2));
        
        outSigOut = obj.processedSignal - baseline';
        outMessage = [ 'Baseline removed, method: ' method];
    case 'poly'
        for iChannel = 1:obj.nChannel
           if isempty(parameter); N=floor(2+obj.nSamp/(20*obj.samp)); else N=parameter;end
            index = 1:obj.nSamp;
            [p] = polyfit(index, obj.processedSignal(iChannel,:),N);
            
            outSig = obj.processedSignal(iChannel,:)-polyval(p,index);
            outMessage =  [ 'Baseline removed, method: ' method ' order: ' num2str(N)];
            outSigOut(iChannel,:) = outSig;
            
        end
    case 'filter'
        if isempty(parameter); parameter=[.5 150];end
        outSigOut = obj.applyBandPassFilter(parameter,400);
        outMessage = [ 'Baseline removed, method: ' method ];
    case 'wavelet'
        if isempty(parameter); parameter=[.5 150];end
        outSigOut = obj.applyWaveletFilter(20,parameter);
        outMessage = [ 'Baseline removed, method: ' method];
        
    case 'spline'
        if isempty(parameter(1)); parameter(1) = 670;end
        margin = floor(parameter(1)*obj.samp/1000); %  No of ms before the R wave where to select the knot.
        
        Sig = obj.removeBaseline('savitzkygolay');
        Sig(obj.badLeads==1,:)=[];
        RMSV = sqrt(mean((Sig.^2)));    
        [x, indexPosition] = fFindPeaks(RMSV,'MINPEAKDISTANCE',parameter(2));      
        refSize = floor(20*obj.samp/1000); % size of the reference interval.
        
        
        nInt = numel(indexPosition);
        refIntervals = cell(1,2);
        %refIntervals{1} = zeros(1,nInt);
        %refIntervals{2} = zeros(1,nInt);
        falseStart = false;
        k=1;
        for r = 1:nInt
            
            tmp = indexPosition(r) - margin;
            if tmp >= 1
                refIntervals{1}(k) = tmp;
                refIntervals{2}(k) = tmp + refSize;
                k = k+1;
            else
                falseStart = true;
            end
            
        end
        %% For each lead
        for l = 1:obj.nChannel
            %% Average intercals to find reference points
            NR = numel(refIntervals{2});
            refPoints{l} = zeros(2,NR);
            for nr = 1:NR
                %nr
                refPoints{l}(1,nr) = mean( obj.processedSignal(l, refIntervals{1}(nr):refIntervals{2}(nr)));
                refPoints{l}(2,nr) = (refIntervals{2}(nr) + refIntervals{1}(nr) )/ 2;
            end
            
            %% fit splines to reference points to obtain baseline
            baseline(l,:) = spline(refPoints{l}(2,:),refPoints{l}(1,:), [1:1:obj.nSamp]);
            
            %% substract baseline from the original signal
            outSigOut(l,:) = obj.processedSignal(l,:) - baseline(l,:);
            
            outMessage = [ 'Baseline removed, method: ' method];

        end%end of all leads
end



% Output
if nargout==0
    obj.pProcessedSignal = outSigOut;
    obj.filteringDetails{size(obj.filteringDetails,2)+1} =  outMessage;
    % Plot
    if (obj.pdisplay==1)
        plot(obj);
    end
else
    if (obj.pdisplay==1)
        plot(outSigOut);
    end
end
