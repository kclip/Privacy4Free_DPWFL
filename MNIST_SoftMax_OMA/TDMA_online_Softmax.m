

% Initialize w
W_TDMA=W0;






gamma=Bound_gamma;
     



alpha_off_TDMA=zeros(T/K,K);

res_dp=ones(1,K)*dp;

for iter=1:T/K
    
      grad_inv= - kron(ones(1, 10 ), X) .* kron( ( Y- exp(X*W_TDMA)./sum(exp(X*W_TDMA),2) ) ,ones(1,d) );
    
      grad_inv= min(1, Bound_gamma./sqrt(sum(grad_inv.^2,2))).* grad_inv;
%      gamma=max(sqrt(sum(grad_inv.^2,2)));
    
   for k=1:K
       
        h_k=h((iter-1)*K+k,k);
        
        %j-i => (t'-t)*K
        %input h^2, from current iter to T/K
        pred_h=(kappa+rho.^(2*[0:T/K-iter]*K))/(kappa+1)*h_k^2+(1-rho.^(2*[0:T/K-iter]*K))/(kappa+1);
       
       grad_k=sum(grad_inv(sum(D(1:k))-D(k)+1:sum(D(1:k)),:),1)'/D(k)+ lambda*2*reshape(W_TDMA,prod(size(W_TDMA)),1);
  
%         if iter==1
%             G_k=gamma;
%         else
%             G_k=norm(recv_grad(:,k))/D(k);
%         end
        G_k=norm(grad_k);
        
        
        
        %%% Optimal PA %%%%
        alpha_min=sqrt(P)/D(k)/G_k;

        if sum( (sqrt(2)*gamma)^2*P*pred_h/(sigma_noise*D(k)*G_k)^2 )<res_dp(k)
            alpha_k=alpha_min;
        else
            lb_sqrt_lambda=0;
            ub_sqrt_lambda=10^8;
  
            while ub_sqrt_lambda-lb_sqrt_lambda >10^(-3)

                mid_lambda_TDMA=(lb_sqrt_lambda+ub_sqrt_lambda)/2;

                if sum( (sqrt(2)*gamma)^2*  min( (1-mu/L).^-([iter:T/K]/2)/(sqrt(2*(1-mu/L)^(-T/K)*mid_lambda_TDMA)*gamma), P*pred_h/(sigma_noise*D(k)*G_k)^2) ) <=res_dp(k)
                    ub_sqrt_lambda=mid_lambda_TDMA;
                else
                    lb_sqrt_lambda=mid_lambda_TDMA;        
                end
            end
             mid_lambda_TDMA=(lb_sqrt_lambda+ub_sqrt_lambda)/2;
            alpha_k=min((1-mu/L)^(-iter/4)*sigma_noise*(2*(1-mu/L)^(-T/K)*mid_lambda_TDMA)^(-1/4)/h_k/sqrt(gamma), alpha_min);
            
            rec_sol_TDMA(iter,T_index,ave_index)=mid_lambda_TDMA;
        end  
        

        %%% Optimal PA end %%%
        
        
%        if iter==1
%             alpha_off_TDMA(:,k)=min((1-mu/L).^(-[iter:K/T]/4)*sigma_noise*(2*mid_lambda)^(-1/4)./sqrt(pred_h*gamma), alpha_min);
%        end
       
       recv_grad(:,k)=D(k)*grad_k+ channel_noise((iter-1)*K+k,:)'/(alpha_k*h_k);    %channel_noise=sigma_noise*randn(T,K,d);
       
       res_dp(k)=res_dp(k)-(sqrt(2)*gamma*h_k*alpha_k)^2/sigma_noise^2;
       
   end
    
   clear grad_inv
   
  
   W_TDMA=W_TDMA-eta*reshape(sum(recv_grad,2)/n,d,10);
   

    %projection w
    norm_W=sqrt(sum(sum(W_TDMA.^2)));
    if norm_W>max_w
        W_TDMA=W_TDMA/norm_W*max_w;
    end
%     Loss_iter_TDMA(iter)= -sum(sum(Y.*log( exp(X*W_TDMA)./sum(exp(X*W_TDMA),2) )))/n+lambda*sum(sum(W_TDMA.^2));
% 
%     recod_grad_norm(iter)=sqrt(sum((sum(recv_grad,2)/n).^2));
% %     recod_w_norm(iter)=sqrt(sum(sum(W.^2)));
%     
%     [M,pred_label]= max(exp(X_test*W_TDMA)./sum(exp(X_test*W_TDMA),2),[],2);
%     pred_label=pred_label-1;
%     pred_err(iter)= sum ( test_set(:,1)~= pred_label) /size(test_set,1);
%  
end

% 
% figure(1)
% plot(Loss_iter_TDMA) 
% figure(2)
% plot(pred_err)
% figure(3)
% plot(recod_grad_norm)
%    
