

% Initialize w
w_TDMA=zeros(d,1);





gamma=0;

alpha_off_TDMA=zeros(T/K,K);

res_dp=ones(1,K)*dp;

for iter=1:T/K
    
     grad_inv=(Y-X*w_TDMA).*(-X);
     grad_inv= diag(min(1, Bound_gamma./sqrt(sum(grad_inv.^2,2))))*grad_inv;
%      gamma=max(sqrt(sum(grad_inv.^2,2)));
     gamma=Bound_gamma;
     
     G_k=0;
    
   for k=1:K
       
        h_k=h((iter-1)*K+k,k);
        
        %j-i => (t'-t)*K
        %input h^2, from current iter to T/K
        pred_h=(kappa+rho.^(2*[0:T/K-iter]*K))/(kappa+1)*h_k^2+(1-rho.^(2*[0:T/K-iter]*K))/(kappa+1);
       
        grad_k=sum(grad_inv(sum(D(1:k))-D(k)+1:sum(D(1:k)),:),1)'/D(k)+ lambda*2*w_TDMA;

        if iter==1
            G_k=gamma;
        else
            G_k=norm(recv_grad(:,k))/D(k);
        end
        
        
        
        %%% Optimal PA %%%%
        alpha_min=sqrt(P)/D(k)/G_k;

        if sum( (sqrt(2)*gamma)^2*P*pred_h/(sigma_noise*D(k)*G_k)^2 )<res_dp(k)
            alpha_k=alpha_min;
        else
            lb_sqrt_lambda=0;
            ub_sqrt_lambda=10^8;
  
            while ub_sqrt_lambda-lb_sqrt_lambda >10^(-8)

                mid_lambda_TDMA=(lb_sqrt_lambda+ub_sqrt_lambda)/2;

                if sum( (sqrt(2)*gamma)^2*  min( (1-mu/L).^-([iter:T/K]/2)/(sqrt(2*(1-mu/L)^(-T/K)*mid_lambda_TDMA)*gamma), P*pred_h/(sigma_noise*D(k)*G_k)^2) ) <=res_dp(k)
                    ub_sqrt_lambda=mid_lambda_TDMA;
                else
                    lb_sqrt_lambda=mid_lambda_TDMA;        
                end
            end
             mid_lambda_TDMA=(lb_sqrt_lambda+ub_sqrt_lambda)/2;
            alpha_k=min((1-mu/L)^(-iter/4)*sigma_noise*(2*(1-mu/L)^(-T/K)*mid_lambda_TDMA)^(-1/4)/h_k/sqrt(gamma), alpha_min);
        end  
        %%% Optimal PA end %%%
        
        
%        if iter==1
%             alpha_off_TDMA(:,k)=min((1-mu/L).^(-[iter:K/T]/4)*sigma_noise*(2*mid_lambda)^(-1/4)./sqrt(pred_h*gamma), alpha_min);
%        end
       
       recv_grad(:,k)=D(k)*grad_k+ channel_noise((iter-1)*K+k,:)'/(alpha_k*h_k);    %channel_noise=sigma_noise*randn(T,K,d);
       
       res_dp(k)=res_dp(k)-(sqrt(2)*gamma*h_k*alpha_k)^2/sigma_noise^2;
       
    end
   
  
   w_TDMA=w_TDMA-mu/L* sum(recv_grad,2)/n;      
   
  %projection w
    if norm(w_TDMA)>max_w
        mark_TDMA=1;
        w_TDMA=w_TDMA/norm(w_TDMA)*max_w;
    end
   
%    Loss_iter_TDMA(iter*K,NumUser_index,ave_index)=( sum((X*w_TDMA-Y).^2)/2/n+5*10^(-5)*norm(w_TDMA)^2 )/L_opt;
  
end



% figure(2)
% plot(Loss_iter_TDMA)
% pause(0.1)
