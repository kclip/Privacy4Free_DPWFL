% Initialize w
W_OA=W0;



G_k=0;

gamma=0;


for iter=1:T/K
    
 
   for k=1:K
       
        h_k=h((iter-1)*K+k,k);
        
                
        index_sel=sum(D(1:k))-D(k)+1:sum(D(1:k));
        grad_k= reshape( - X(index_sel,:)'*( Y(index_sel,:)- exp(X(index_sel,:)*W_OA)./sum(exp(X(index_sel,:)*W_OA),2) ) /D(k) +lambda*2*W_OA , prod(size(W_OA)),1); 

        G_k=norm(grad_k);
        
        
        %%%% Static PA %%%%
        alpha_k=sqrt(P)/D(k)/G_k;

       
        recv_grad(:,k)=D(k)*grad_k+ channel_noise((iter-1)*K+k,:)'/(alpha_k*h_k);    %channel_noise=sigma_noise*randn(T,K,d);
       
       
   end
   clear grad_k

   W_OA=W_OA-eta* reshape(sum(recv_grad,2)/n,d,10);

   %projection w
   norm_W=sqrt(sum(sum(W_OA.^2)));
    if norm_W>max_w
        W_OA=W_OA/norm_W*max_w;
    end
   
% 
%      Loss_iter_OA(iter)= -sum(sum(Y.*log( exp(X*W_OA)./sum(exp(X*W_OA),2) )))/n+lambda*sum(sum(W_OA.^2));
% 
%     recod_grad_norm(iter)=sqrt(sum((sum(recv_grad,2)/n).^2));
% %     recod_w_norm(iter)=sqrt(sum(sum(W.^2)));
%     
%     [M,pred_label]= max(exp(X_test*W_OA)./sum(exp(X_test*W_OA),2),[],2);
%     pred_label=pred_label-1;
%     pred_err(iter)= sum ( test_set(:,1)~= pred_label) /size(test_set,1);
end


% 
% figure(1)
% plot(Loss_iter_OA) 
% figure(2)
% plot(pred_err)
% figure(3)
% plot(recod_grad_norm)