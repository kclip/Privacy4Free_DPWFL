w_NA=zeros(d,1);



G_k=0;

gamma=0;

record_w_NA=[];



for iter=1:T
    grad_inv=(Y-X*w_NA).*(-X);
    gamma=max(sqrt(sum(grad_inv.^2,2)));
  
    h_current=h(iter,:);
       
       
    grad_k=zeros(d,K);
    Scal_grad_k=zeros(d,K);
    DkG_k=0;
    for k=1:K
        grad_k(:,k)=sum(grad_inv(sum(D(1:k))-D(k)+1:sum(D(1:k)),:),1)'/D(k)+ lambda*2*w_NA;
        DkG_k(k)=D(k)*norm(grad_k(:,k));
        Scal_grad_k(:,k)=grad_k(:,k)*D(k);
    end
        
    c_NA= sqrt(P)*min(h_current./DkG_k);
    
    recv_grad=(sum(Scal_grad_k,2)+channel_noise(iter,:)'/c_NA)/n;
    
    w_NA=w_NA-mu/L*recv_grad;
    
    %projection w
    if norm(w_NA)>max_w
        mark_TDMA=1;
        w_NA=w_NA/norm(w_NA)*max_w;
    end
    
%     record_w_NA=[record_w_NA,w_NA];
    
    
%      Loss_iter_NA(iter,NumUser_index,ave_index)=( sum((X*w_NA-Y).^2)/2/n+5*10^(-5)*norm(w_NA)^2 )/L_opt;
    
end  

% figure(6)
% plot(Loss_iter_OA)
% pause(0.1)