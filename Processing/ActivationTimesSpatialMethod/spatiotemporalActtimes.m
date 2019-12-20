%% HELP:
%
%	This function computes the activation times in a cardiac sequence using
%	a spatiotemporal approach.
%	In this spatiotemporal approach, the objective is to find the minimum
%	dVdt weighted by the norm of the spatial gradient at each time
%	instance.
%
%	INPUT:
%		- X - <N,T>double - sequence of potentials on the heart.
%		- D - <3N,N>double - estimator of spatial gradient of potentials.
%		It should be organized such that the first N columns correspond to
%		the gradient in X direction of each node, the following N columns
%		are Y and the last Z.
%		- window - int - (OPTIONAL) window of time over which the temporal
%		derivative is calculated.
%		- alpha - double [0,1] - weight to be used for the gradient. If 1
%		it is direct multiplication, using 0 implies no weighting.
%
%


function [tau,objfun,DX,dXdt,timeconfidence] = spatiotemporalActtimes(X,D,varargin)

	%% parse inputs
		if(size(varargin,2)>0),window=varargin{1};else window=1;end;
% 		if(size(varargin,2)>1),alpha=varargin{2};else alpha=1;end;
        if(size(varargin,2)>1),smootheddiffed=varargin{2}; window=0; end;

		[N,T] = size(X);

	%% Compute gradient
		SS = [eye(N), eye(N), eye(N)];		% add up X,Y,Z entries of gradient
		
		DX = sqrt( SS*(D(:,1:N)*X).^2 );		% compute gradient norm
		DX = DX./repmat(max(DX,[],2),[1,T]);	% normalize to 1 (over time)
		DX(isnan(DX)) = 0;						% eliminate outliers

		
	%% Compute min dVdT (Jaume-approach) if not given already
        if window==0;
            dXdt = smootheddiffed;
        else
            [~, dXdt] = findMinDVDT( X, window, 2);
        end
	
	%% Weight dVdt with spatial gradient norm
		objfun= DX.*dXdt; % objfun now contains the gradient norms weighted by negative temporal slope

	%% FIND MINIMUM OF THE OBJECTIVE FUNCTION
		[val,tau]=min(objfun,[],2);
        
    %% Visualization for debug
    if false
        for i=1:N
            plot(X(i,:)./max(abs(X(i,:))),'r'); hold on;
            plot(dXdt(i,:),'b');
            plot(objfun(i,:),'k');
            legend({'Signal','dXdt','objfun'});
            hold off;
            pause;
        end
    end
        
    %% Addition Matthijs (2016-11-21) for time confidence according to duchateau2016
    timeconfidence = nan(N,1);
    for i = 1:N;
        curpotsdiffed = dXdt(i,:);
        belowzeroindices = find(curpotsdiffed<0);
        auc = sum(curpotsdiffed(belowzeroindices));
        slopesurrogate = dXdt(i,tau(i)); % just the temporal slope <-- use this one, as AUC is also just temporal (for now)
    %     slopesurrogate = val(i); % the 'spatiotemp' slope
        timeconfidence(i) = slopesurrogate/auc;
    %     plot(X(i,:)); hold on; scatter(belowzeroindices,zeros(size(belowzeroindices))); plot(dXdt(i,:),'r');
    %     title(['Slope: ' num2str(slopesurrogate) ' / auc: ' num2str(auc) ' = score: ' num2str(timeconfidence(i))]); hold off;
    %     pause;
    end


end
