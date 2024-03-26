function [X_den,P1,P2,iter] = denoise_TV_MT(Xobs, lambda, P1_init, P2_init, pars)
%This function implements the FISTA method for JTV denoising problems. 
%
% INPUT
% Xobs ..............................observed noisy images. [H, W, T]
% lambda ........................ parameter
% pars.................................parameters structure
% pars.MAXITER ..................... maximum number of iterations
%                                                      (Default=100)
% pars.epsilon ..................... tolerance for relative error used in
%                                                       the stopping criteria (Default=1e-4)
% pars.print ..........................  1 if a report on the iterations is
%                                                       given, 0 if the  report is silenced
% pars.tv .................................. type of total variation
%                                                      penatly.  'iso' for isotropic (default)
%                                                      and 'l1' for nonisotropic
%  
% OUTPUT
% X_den ........................... The solution of the problem 
%                                            min{||X-Xobs||^2+2*lambda*TV(X)}

%Define the Projection onto the box
% if((l==-Inf)&&(u==Inf))
%     project=@(x)x;
% elseif (isfinite(l)&&(u==Inf))
%     project=@(x)(((l<x).*x)+(l*(x<=l)));
% elseif (isfinite(u)&&(l==-Inf))
%      project=@(x)(((x<u).*x)+((x>=u)*u));
% elseif ((isfinite(u)&&isfinite(l))&&(l<u))
%     project=@(x)(((l<x)&(x<u)).*x)+((x>=u)*u)+(l*(x<=l));
% else
%     error('lower and upper bound l,u should satisfy l<u');
% end

% Assigning parameres according to pars and/or default values
flag=exist('pars', 'var');
if (flag&&isfield(pars,'MAXITER'))
    MAXITER=pars.MAXITER;
else
    MAXITER=100;
end
% if (flag&&isfield(pars,'epsilon'))
%     epsilon=pars.epsilon;
% else
%     epsilon=1e-4;
% end
% if(flag&&isfield(pars,'print'))
%     prnt=pars.print;
% else
%     prnt=1;
% end
if(flag&&isfield(pars,'type'))
    type=pars.type;
else
    type='iso';
end

[m,n,T]=size(Xobs);
% clear P; clear R;
if(isempty(P1_init))
    P1 = zeros(m-1,n,T); P2 = zeros(m,n-1,T);
    R1 = zeros(m-1,n,T); R2 = zeros(m,n-1,T);
else
    P1=P1_init;    P2=P2_init;
    R1=P1_init;    R2=P2_init;

end
% tk=1;
tkp1 = 1;
% count=0;
i = 0;
% fval=inf;
% fun_all=[];
% epsilon=1e-4;
D=zeros(m,n,T);%fval=inf;fun_all=[];

switch type
    case 'iso'
        proj_func = @ (p1, p2) proj_iso(p1, p2, m, n);
    case 'l1'
        proj_func = @ (p1, p2) proj_l1(p1, p2);
    otherwise
        error('unknown type of total variation. should be iso or l1');
end

%while((i<MAXITER)&&(count<5))
while(i<MAXITER)    
    %fold=fval;  
    %%%%%%%%%
    % updating the iteration counter
    i = i+1;
    %Dold=D;    
    P1_old = P1;P2_old = P2;    
    tk = tkp1;
%     W = pars.w;   
    %%%%%%%%%%
    %Computing the gradient of the objective function
    D = Xobs-lambda*Lforward(R1, R2, m, n);

    %%%%%%%%%%
    % Taking a step towards minus of the gradient
    step0=1/(8*lambda);
	[G1, G2] = Ltrans(D, m, n);
    P1 = R1 + step0.*G1;
    P2 = R2 + step0.*G2;
    %%%%%%%%%%
    % Peforming the projection step
    [P1, P2] = proj_func(P1, P2);
    %%%%%%%%%%
    %Updating R and t
    tkp1 = (1+sqrt(1+4*tk^2))/2;
    
    step = (tk-1)/tkp1;

    R1 = P1 + step.*(P1-P1_old);
    R2 = P2 + step.*(P2-P2_old);
        
%     dD = D-Dold;
%     norm_dD = norm(dD(:),'fro');
%     norm_D = norm(D(:),'fro');
%     re=norm_dD/norm_D;
% 
%     if (re<epsilon)
%         count=count+1;
%     else
%         count=0;
%     end
% 
%     C = Xobs-lambda*Lforward(P1, P2, m, n);
%     fval=norm(C(:),'fro')^2;
% %    
%     fun_all=[fun_all;fval];
%     fprintf('          %5d                   %10.10f              %10.10f',i,fval,re);
%     if (fval>fold)
%         fprintf('  *\n');
%     else
%         fprintf('   \n');
%     end    
end
X_den = D;
iter = i;

function [P1, P2] = proj_iso(P1, P2, m, n)
    A = zeros(m,n,size(P1,3));
    A(1:m-1,:,:) = P1.^2;
    A(:,1:n-1,:) = A(:,1:n-1,:) + P2.^2;
    A=sqrt(max(A,1));
   
    P1 = P1./A(1:m-1,:,:); 
    P2 = P2./A(:,1:n-1,:);
    
end

function [P1, P2] = proj_l1(P1, P2)
    P1=P1./(max(abs(P1),1));
    P2=P2./(max(abs(P2),1));
end
        
end