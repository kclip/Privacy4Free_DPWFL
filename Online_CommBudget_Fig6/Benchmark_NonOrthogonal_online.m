w_NAbm=zeros(d,1);



G_k=0;

gamma=Bound_gamma/5;



for iter=1:T
    grad_inv=(Y-X*w_NAbm).*(-X);
    grad_inv= diag(min(1, gamma./sqrt(sum(grad_inv.^2,2))))*grad_inv;

%     gamma=max(sqrt(sum(grad_inv.^2,2)));
    
    
    h_current=h(iter,:);
       
       
    grad_k=zeros(d,K);
    Scal_grad_k=zeros(d,K);
    
    if iter==1 
        DkG_k=ones(1,K)*n/K*gamma;
    else
        DkG_k= ones(1,K)*norm(recv_grad*n/K);
    end
        
    
    c_NAbm=min(sqrt(sigma_noise^2*dp/2/T/gamma^2), sqrt(P)*min(h_current./DkG_k));
    
    
     for k=1:K
        grad_k(:,k)=sum(grad_inv(sum(D(1:k))-D(k)+1:sum(D(1:k)),:),1)'/D(k)+ lambda*2*w_NAbm;
%         DkG_k(k)=D(k)*norm(grad_k(:,k));

        Scal_grad_k(:,k)=grad_k(:,k)*D(k)*min(1, sqrt(P)/(c_NAbm/h_current(k))/norm(grad_k(:,k)*D(k)));
    end
    
    recv_grad=(sum(Scal_grad_k,2)+channel_noise(iter,:)'/c_NAbm)/n;
    
    w_NAbm=w_NAbm-mu/L*recv_grad;
    
    %projection w
    if norm(w_NAbm)>max_w
        mark_NAbm=1;
        w_NAbm=w_NAbm/norm(w_NAbm)*max_w;
    end
    
%     Loss_iter_NAbm(iter,NumUser_index,ave_index)=( sum((X*w_NAbm-Y).^2)/2/n+5*10^(-5)*norm(w_NAbm)^2 )/L_opt;
    
end  

% figure(6)
% plot(Loss_iter_OAbm)
% pause(0.1)

   