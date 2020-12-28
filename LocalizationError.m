function [LE, invPacing] = LocalizationError(ATrec,TriRec,ATinv,TriInv,pacingSite)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    %  for i = 1:size(TriInv.X,1) % Apply median filter to ATs with neighbours
    %      [neighbour, distance] = findNearNeighbours(TriInv,i);
    %      ATinvAvg(i) = median(ATinv([i neighbour']));
    %  end
    %
    
if(isempty(pacingSite))
    [~, pacingSite] = min(ATrec);
end


%% Define if single point or group of points have earliest activation
[ATinvOrdered Order] = sort(ATinv);
[earlySites] = find(ATinvOrdered==ATinvOrdered(1));
if(length(earlySites)==1) %  Single early activation
    invPacing =  TriInv.X(Order(1),:);
else % Multiple points with same activation, try to group them
    group{1} = Order(1);
    k=1;n=1;
    Order = Order((1:length(earlySites)));
    OrderAll = Order;
    while(~isempty(Order))
        [neighbour, distance] = findNearNeighbours(TriInv,Order(1));
        neighbour = intersect(neighbour,OrderAll);
        if(isempty(neighbour))
            % single point with no neighbours
            k = k+1;
            group{k} =  Order(1);
            Order(1)=[];
        else
            if(isempty(intersect([neighbour; Order(1)],group{k})))
                % check other groups also not connected
                m = k-1;
                while(m>0)
                    %disp('second while loop')
                    if(isempty(intersect([neighbour; Order(1)],group{m})))
                        m = m-1;
                    else
                        group{m} = unique([group{m} neighbour']);
                        Order(1)=[];
                        m=-1;
                    end
                end
                if(m==0)
                    k = k+1;
                    group{k} =  Order(1);
                    Order(1)=[];
                end
            else
                group{k} = unique([group{k} neighbour']);
                Order(1)=[];
            end
        end
    end
    % Defines pacing site as the centre of the biggest group
    if(size(group,2)==1)
        if(size(group{1},2)>1)
            invPacing =  mean(TriInv.X(group{1},:));
        else
            invPacing =  TriInv.X(group{1},:);
        end
    else % multiple groups; define biggest
        [~,id] = max(cellfun(@numel,group));
        if(size(group{id},2)>1)
            invPacing =  mean(TriInv.X(group{id},:));
        else
            invPacing =  TriInv.X(group{id},:);
        end
        %invPacing =  mean(TriInv.X(group{id},:));
    end
        
end


    
%% Simple Localization Error
%invPacing
LE = EuclideanDistance(invPacing,TriRec.X(pacingSite,:));
 
%% output AT maps
% Dessin(TriInv,ATinv,ATinv)
% hold on
% %invPacing
% scatter3(invPacing(1),invPacing(2),invPacing(3),'ko')   
%     
return

