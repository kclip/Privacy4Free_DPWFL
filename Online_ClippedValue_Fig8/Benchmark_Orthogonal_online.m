

% Initialize w
w_OAbm=zeros(d,1);



G_k=0;

gamma=0;


for iter=1:T/K
    
     grad_inv=(Y-X*w_OAbm).*(-X);
     grad_inv= diag(min(1, Bound_gamma./sqrt(sum(grad_inv.^2,2))))*grad_inv;
     
%      gamma=max(sqrt(sum(grad_inv.^2,2)));
     gamma=Bound_gamma;
     G_k=0;

   for k=1:K
       
        h_k=h((iter-1)*K+k,k);
        
        grad_k=sum(grad_inv(sum(D(1:k))-D(k)+1:sum(D(1:k)),:),1)'/D(k)+ lambda*2*w_OAbm;

        if iter==1
            G_k=gamma;
        else
            G_k=norm(recv_grad(:,k))/D(k);
        end
        
        %%%% Static PA %%%%
        alpha_k=min(sqrt(sigma_noise^2*dp/2/T/(h_k*gamma)^2),sqrt(P)/D(k)/G_k);

       
        recv_grad(:,k)=D(k)*grad_k+ channel_noise((iter-1)*K+k,:)'/(alpha_k*h_k);    %channel_noise=sigma_noise*randn(T,K,d);
       
       
   end
   
   w_OAbm=w_OAbm-mu/L* sum(recv_grad,2)/n;
   
   %projection w
    if norm(w_OAbm)>max_w
        mark_OAbm=1;
        w_OAbm=w_OAbm/norm(w_OAbm)*max_w;
    end
   
%     Loss_iter_OAbm(iter*K,NumUser_index,ave_index)=( sum((X*w_OAbm-Y).^2)/2/n+5*10^(-5)*norm(w_OAbm)^2 )/L_opt;
  
end

% 
% figure(5)
% plot(Loss_iter_OAbm)
