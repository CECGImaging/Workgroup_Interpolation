function smoothedacttimes=smoothactivationtimes(L,activationtimes,lambdas,varargin)


Y=activationtimes;
A=eye(size(L,2));

MANUALFLAG=0;
if(numel(varargin)>0)
    if(strcmpi(varargin{1},'manual'))
        MANUALFLAG=1;
    end
end

LL=(L'*L);
AA=(A'*A);

if length(lambdas)==1
    doManualLambda = true; 
else
    doManualLambda = false;
end

for i=1:size(Y,2)
    
    if doManualLambda
        lambda = lambdas;
    else
        eta=zeros(numel(lambdas),1);
        rho=zeros(numel(lambdas),1);
        for k=1:numel(lambdas)
            temp = (AA+lambdas(k)*LL)\(A'*Y(:,i));
            eta(k)=norm(L*temp);
            rho(k)=norm(Y(:,i)-A*temp);
        end

        [lambda,rho_c,eta_c] = l_corner(rho,eta,lambdas);
        figure; loglog(rho,eta); hold on;  loglog(rho_c,eta_c,'ro'),hold off
    %     
    %     if(MANUALFLAG==1)
    %         figure(1), loglog(rho,eta)
    %         [rho_c,eta_c]=ginput(1); [val,ind]=min(sum([rho-rho_c,eta-eta_c].^2,2)); lambda=lambdas(ind);
    %         hold on, loglog(rho_c,eta_c,'ro'),hold off
    %     else
    %         rholog=log10(rho); etalog=log10(eta); Trho=2*numel(lambdas);
    %         etalog=spline(rholog,etalog,linspace(min(rholog),max(rholog),Trho));
    %         etalog=lowpassma(etalog,10);
    %         lambdaspline=spline(rholog,lambdas,linspace(min(rholog),max(rholog),Trho));
    %         rholog=linspace(min(rholog),max(rholog),Trho);
    %         detalog=diff(etalog,1);
    %         [~,minslopeind]=min(detalog);
    % 
    %         lambda=lambdaspline(minslopeind);
    %     end
    %     
    %     figure(1)
    %     loglog(rho,eta,'b')
    %     hold on
    %     loglog(10^rholog(minslopeind),10^etalog(minslopeind),'ro')
    %     hold off
    %     pause
    end
    
    X_reg(:,i) = (AA+lambda*LL)\(A'*Y(:,i));
    
    fprintf('t=%i\n',i)
    
    disp(['Used lambda: ' num2str(lambda)]);
    
end

smoothedacttimes=X_reg;  
