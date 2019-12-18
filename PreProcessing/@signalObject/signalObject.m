classdef signalObject < handle
    % Custom class to manage signals 
    %   Signal Object provides methods to manipulate a matrix of signals.
    %   
    % Input: 
    % - signalMatrix with the data, rows are the channels, and column
    % the time samples
    % 
    %Optional inputs:
    % - 'label': label 
    % - 'display': boolean value, to plot/or not signals when processing is perfomed
    % - 'tri': TriRep geometry associated to the signals
    % - 'signalType': 'unipolar' or 'bipolar' usefull to apply different processing according to the signal type (fft for example)
    % - 'samp': sampling rate in HZ (default 1kHz)
    % - 'badLeads': boolean vector with the bad leads
    % 
    
    
    
    properties (Dependent)
        rawSignal;  % Signal
        processedSignal; % Processed Signal
        tri = '';    % TriRep if any

        
    end
    properties (GetAccess = 'public', SetAccess = 'public')
        badLeads;        % if =0, good lead
        sigLabel='';   % channel label
        markers=[];
        markersLabel=[];
        triRepIndex = 0; % assign X index in triRep function
        signalType = 'unipolar';   % Signal type
    end
    
    properties (GetAccess = 'public', SetAccess = 'private')
%        tri = '';    % TriRep if any
        nChannel;
        nSamp; % Number of samples
        samp;   % Sampling rate
        filteringDetails = ''; % string that contains the filtering details
    end
    
    properties (Access = 'private')
        pTri;
        pRawSignal;                                   % Raw signal
        pProcessedSignal;            % filteredSignal
        pdisplay = 0;                       % Display options
        
    end
    
    
    
    methods
        %% Constructor
        function obj = signalObject(signal,varargin)
            
            obj.pRawSignal = signal;
            obj.pProcessedSignal = signal;
            obj.nSamp = size(signal,2);
            obj.nChannel = size(signal,1);
            
            opts = parseArgsToStruct(varargin, ...
                {'label', 'display','tri','signalType','samp','badLeads'}, ...
                {'sigLabel', 'pdisplay','tri','signalType','samp','badLeads'}, ...
                {'',  0,[], 'unipolar',1000,0 });
            
            fnames = fieldnames(opts);
            for i = 1:length(fnames), obj.(fnames{i}) = opts.(fnames{i}); end
            if obj.badLeads==0,             obj.badLeads = false(1,size(signal,1)); end
            
        end
        
        
        
        %% Get/Set
        function outraw = get.tri(obj)
            outraw = obj.pTri;
        end
  
        function outraw = get.processedSignal(obj)
            outraw = obj.pProcessedSignal;
        end
        function resetProcessedSignal(obj)
            % reset all the processing perfomed and restored the raw data
            obj.pProcessedSignal = obj.pRawSignal;
            obj.filteringDetails = '';
        end
        function outraw = get.rawSignal(obj)
            outraw = obj.pRawSignal;
        end
        
        function set.tri(obj,tri)
            if ~isempty(tri)
                if isa(tri,'TriRep')
                    obj.pTri = tri;
                else
                    error('tri coordinates not valides');
                end
            end
        end
        
        function set.signalType(obj,sigType)
            switch lower(sigType)
                case {'unipolar','bipolar','ap'}
                    obj.signalType=lower(sigType);
                otherwise
                    error('Value must be either unipolar, bipolar or AP');
            end
        end
        %%
        function minusConstant(obj,value,label)
            % Removes a constant vector or a constant value to all the channels
            if nargin==2, label=''; else, label = [label ' '];end
            if((size(value,1) == 1)  && (size(value,2) == size(obj.processedSignal,2)) )
                obj.pProcessedSignal = obj.processedSignal-repmat(value,obj.nChannel,1);
            elseif((size(value,1) == 1)  && (size(value,2) == size(obj.processedSignal,1)) )
                obj.pProcessedSignal = obj.processedSignal-repmat(value,size(obj.processedSignal,2),1)';
            elseif((size(value,1) == size(obj.processedSignal,1))  && (size(value,2) == size(obj.processedSignal,2)) )
                obj.pProcessedSignal = obj.processedSignal-value;
            end
            obj.filteringDetails{size(obj.filteringDetails,2)+1} =  [ label 'value Removed' ];
            
        end
        
         %%
        function minusAvg(obj,indexChan,label)
            % Removes a average of the indexChan to all the signals
            if nargin<=2, label=''; else, label = [label ' '];end
            if nargin<=1, indexChan = find(obj.badLeads==0);end
            avgSignal = mean(obj.processedSignal(indexChan,:));
            if nargin<=1,
                obj.minusConstant(avgSignal,[label 'all']);
            else
                obj.minusConstant(avgSignal,[label sprintf('%d',indexChan)]);
            end
        end
        
        %% add/remove Marker
        function addMarker(obj,position, label)
            % Add temporal marker
            % input:
            % position : [onset offset] (offset optional)
            % label: maker label
            if ~iscell(position), position{1}=position;end
            if ~iscell(label), label{1}=label;end
            
            nMarker = size(obj.markers,1);
            if size(position{1},1)~=1 || size(position{1},2)==0  || size(position{1},2)>2
                error('Need an Onset Marker, an Offset Marker (optional) and a label');
            end
            obj.markers{nMarker+1} = position{1};
            obj.markersLabel{nMarker+1} = label{1};
        end
        
        function removeMarker(obj,index)
            % remove a marker from the list
            % input (optional): index of the marker to remove
            nMarker = size(obj.markers,2);
            if (nargin==1 || isempty(index))
                for iMarker = 1:nMarker
                    disp([ '[' num2str(iMarker)  '] - ' obj.markersLabel{iMarker} ' - [ ' num2str(obj.markers{iMarker}) ']'   ]);
                end
            else
                indexBool = true(1,nMarker);
                indexBool(index) = false;
                obj.markers=obj.markers(indexBool,:);
                cp=1;
                for ind=1:length(indexBool)
                    if indexBool(ind)
                        A{cp}=obj.markersLabel{ind};
                        cp=cp+1;
                    end
                end
                obj.markersLabel = A;
            end
        end
        
        function viewMarker(obj)
            % list the markers
            nMarker = size(obj.markers,2);
            for iMarker = 1:nMarker
                disp([ '[' num2str(iMarker)  '] - ' obj.markersLabel{iMarker} ' - [ ' num2str(obj.markers{iMarker}) ']'   ]);
            end
        end
        
        %% Ploting function
        function plot(obj,channel,indexMarker)
            % Plot function is used to plot the signal either one by one,
            % or all surimposed
            % input:
            % obj: signal Object
            % channel: channel to display if  ommitted or empty
            % displays all channels
            % indexMarker: if single value; it is read as a marker index
            % if more than 1 value, it is read as the temporal index to be
            % plotted
            % 
            % keypress:
            % uparrow/downarrow to navigate in the channels
            % backspace to mark the current channel as badLeads
            % space to mark the current channel as good leads
            
            
            if nargin==1|| isempty(channel) channel = 1:obj.nChannel;end
            
            if nargin==3 && size(indexMarker,2)>=2
                    index = indexMarker;
            elseif nargin==3,
                if size(obj.markers{indexMarker},2)==2
                    index = obj.markers{indexMarker}(1):obj.markers{indexMarker}(2);
                elseif size(obj.markers{indexMarker},2)==1
                    index = obj.markers{indexMarker}(1):obj.markers{indexMarker}(1)+min(500+obj.markers{indexMarker}(1),obj.nSamp);
                end
            else
                index = 1:obj.nSamp;
            end

            % Figure plot
            
            scaleBar = median(max(obj.processedSignal')-min(obj.processedSignal'))/10;%1000;%1/10*std(obj.pProcessedSignal(:));
            plotSig1 = obj.pRawSignal(:,index);
            plotSig2 = obj.pProcessedSignal(:,index);
            if size(channel,2)==1
                hplot1 = []; hplot2 = [];
                hfig = figure('Name', 'Check Filtering - press up/down keys to change channel', ...
                    'KeypressFcn', {@fkeypress}, 'Units', 'Normalized', 'Position', [0.2 0.3 0.6 0.4]);
                hold on
                if(obj.badLeads(channel)==false)
                    hplot1 = plot(plotSig1(channel,:)','b');
                    hplot2 = plot(plotSig2(channel,:)','color',[0 0.2 0]);
                    hScaleBar = line([0 0],[mean(plotSig2(channel(1),:))-scaleBar mean(plotSig2(channel(1),:))+scaleBar],'LineWidth',3,'Color','m');

                    legend('Original','Processed');
                else
                    hplot1 = plot(plotSig1(channel,:)','r');
                    hplot2 = plot(plotSig2(channel,:)','r');
                    hScaleBar = line([0 0],[mean(plotSig2(channel(1),:))-scaleBar mean(plotSig2(channel(1),:))+scaleBar],'LineWidth',3,'Color','m');

                    legend('Bad Channel');
                end
                title(['Channel #' num2str(channel)] );
            else
                plot(plotSig1(channel,index)','b');
                hold on
                plot(plotSig1(channel,index)','r');
                title([obj.sigLabel] );
            end
            
            hold off
            % h = warndlg('Scroll through signals using arrows, press backspace to mark as a bad channel (red) and return to mark as a good channel (blue)','Instructions');
            % waitfor(hfig);
            
            % Key presses
            function fkeypress(hObj,event)
                if strcmpi(event.Key, 'downarrow') || strcmpi(event.Key, 'rightarrow')
                    channel = channel+1;
                    if(channel ==(size(plotSig1,1)+1)); channel =1; end
                    if(obj.badLeads(channel)==false)
                        set(hplot1, 'YData', plotSig1(channel,:)','Color','b');
                        set(hplot2, 'YData', plotSig2(channel,:)','color',[0 0.2 0]);
                        set(hScaleBar,'YData',[mean(plotSig2(channel(1),:))-scaleBar mean(plotSig2(channel(1),:))+scaleBar]);

                        legend('Original','Processed');
                    else
                        set(hplot1, 'YData', plotSig1(channel,:)','Color','r');
                        set(hplot2, 'YData', plotSig2(channel,:)','Color','r');
                         set(hScaleBar,'YData',[mean(plotSig2(channel(1),:))-scaleBar mean(plotSig2(channel(1),:))+scaleBar]);
                       legend('Bad Channel');
                    end
                    title(['Channel #' num2str(channel)] );
                elseif (strcmpi(event.Key, 'uparrow') || strcmpi(event.Key, 'leftarrow'))
                    channel = channel-1;
                    if(channel==0); channel = size(plotSig1,1);end
                    if(obj.badLeads(channel)==false)
                        set(hplot1, 'YData', plotSig1(channel,:)','Color','b');
                        set(hplot2, 'YData', plotSig2(channel,:)','color',[0 0.2 0]);
                        set(hScaleBar,'YData',[mean(plotSig2(channel(1),:))-scaleBar mean(plotSig2(channel(1),:))+scaleBar]);
                        legend('Original','Processed');
                    else
                        set(hplot1, 'YData', plotSig1(channel,:)','Color','r');
                        set(hplot2, 'YData', plotSig2(channel,:)','Color','r');
                        set(hScaleBar,'YData',[mean(plotSig2(channel(1),:))-scaleBar mean(plotSig2(channel(1),:))+scaleBar]);
                        legend('Bad Channel');
                    end
                    title(['Channel #' num2str(channel)] );
                elseif strcmpi(event.Key, 'backspace') || strcmpi(event.Key, 'b')
                    obj.badLeads(channel) = true;
                    set(hplot1, 'YData', plotSig1(channel,:)','Color','r');
                    set(hplot2, 'YData', plotSig2(channel,:)','Color','r');
                    set(hScaleBar,'YData',[mean(plotSig2(channel(1),:))-scaleBar mean(plotSig2(channel(1),:))+scaleBar]);
                    legend('Bad Channel');
                elseif strcmpi(event.Key, 'return') || strcmpi(event.Key, 'space')
                    obj.badLeads(channel) = false;
                    set(hplot1, 'YData', plotSig1(channel,:)','Color','b');
                    set(hplot2, 'YData', plotSig2(channel,:)','color',[0 0.2 0]);
                    set(hScaleBar,'YData',[mean(plotSig2(channel(1),:))-scaleBar mean(plotSig2(channel(1),:))+scaleBar]);
                    legend('Original','Processed');
                end
                
            end
            
            
        end
        
        
    end
end

