

% Initialize w
W_OAbm=W0;


G_k=0;

gamma=Bound_gamma;


for iter=1:T/K
    
    grad_inv= - kron(ones(1, 10 ), X) .* kron( ( Y- exp(X*W_OAbm)./sum(exp(X*W_OAbm),2) ) ,ones(1,d) );
    
    grad_inv= min(1, Bound_gamma./sqrt(sum(grad_inv.^2,2))).* grad_inv;
     


   for k=1:K
       
        h_k=h((iter-1)*K+k,k);
        
        grad_k=sum(grad_inv(sum(D(1:k))-D(k)+1:sum(D(1:k)),:),1)'/D(k)+ lambda*2*reshape(W_OAbm,prod(size(W_OAbm)),1);
        

%         if iter==1
%             G_k=gamma;
%         else
%             G_k=norm(recv_grad(:,k))/D(k);
%         end
%         

        G_k=norm(grad_k);
        %%%% Static PA %%%%
        alpha_k=min(sqrt(sigma_noise^2*dp/2/T/(h_k*gamma)^2),sqrt(P)/D(k)/G_k);

       
        recv_grad(:,k)=D(k)*grad_k+ channel_noise((iter-1)*K+k,:)'/(alpha_k*h_k);    %channel_noise=sigma_noise*randn(T,K,d);
       
       
   end
   
   clear grad_inv

   
   W_OAbm=W_OAbm-eta*reshape(sum(recv_grad,2)/n,d,10);
   
    %projection w
    norm_W=sqrt(sum(sum(W_OAbm.^2)));
    if norm_W>max_w
        W_OAbm=W_OAbm/norm_W*max_w;
    end
%     Loss_iter_OAbm(iter)= -sum(sum(Y.*log( exp(X*W_OAbm)./sum(exp(X*W_OAbm),2) )))/n+lambda*sum(sum(W_OAbm.^2));
% 
%     recod_grad_norm(iter)=sqrt(sum((sum(recv_grad,2)/n).^2));
% %     recod_w_norm(iter)=sqrt(sum(sum(W.^2)));
%     
%     [M,pred_label]= max(exp(X_test*W_OAbm)./sum(exp(X_test*W_OAbm),2),[],2);
%     pred_label=pred_label-1;
%     pred_err(iter)= sum ( test_set(:,1)~= pred_label) /size(test_set,1);
%   
end

% figure(1)
% plot(Loss_iter_OAbm) 
% figure(2)
% plot(pred_err)
% figure(3)
% plot(recod_grad_norm)
%    
