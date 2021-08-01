% Initialize w
w_OACLP=zeros(d,1);



G_k=0;



for iter=1:T/K
    
     grad_inv=(Y-X*w_OACLP).*(-X);
     grad_inv= diag(min(1, Bound_gamma./sqrt(sum(grad_inv.^2,2))))*grad_inv;

     

   for k=1:K
       
        h_k=h((iter-1)*K+k,k);
        
        grad_k=sum(grad_inv(sum(D(1:k))-D(k)+1:sum(D(1:k)),:),1)'/D(k)+ lambda*2*w_OACLP;

        G_k=norm(grad_k);
        
        
        %%%% Static PA %%%%
        alpha_k=sqrt(P)/D(k)/G_k;

       
        recv_grad(:,k)=D(k)*grad_k+ channel_noise((iter-1)*K+k,:)'/(alpha_k*h_k);    %channel_noise=sigma_noise*randn(T,K,d);
       
       
   end
   
   w_OACLP=w_OACLP-mu/L* sum(recv_grad,2)/n;
   
     %projection w
    if norm(w_OACLP)>max_w
        mark_TDMA=1;
        w_OACLP=w_OACLP/norm(w_OACLP)*max_w;
    end
   
%    Loss_iter_OA(iter*K,NumUser_index,ave_index)=( sum((X*w_OA-Y).^2)/2/n+5*10^(-5)*norm(w_OA)^2 )/L_opt;
  
end

% 
% figure(5)
% plot(Loss_iter_OA)