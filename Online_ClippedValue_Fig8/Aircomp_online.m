w_aircomp=zeros(d,1);



G_k=0;

gamma=0;

record_w_aircomp=[];

res_dp=dp;



for iter=1:T
    grad_inv=(Y-X*w_aircomp).*(-X);
    
    grad_inv= diag(min(1, Bound_gamma./sqrt(sum(grad_inv.^2,2))))*grad_inv;


%     gamma=max(sqrt(sum(grad_inv.^2,2)));
    
    gamma=Bound_gamma;
    
    h_current=h(iter,:);
        
       
    %input h^2, from current iter to T
    pred_h=kron((kappa+rho.^(2*[0:T-iter]'))/(kappa+1),h_current.^2)+kron((1-rho.^(2*[0:T-iter]'))/(kappa+1),ones(1,K));
       
    grad_k=zeros(d,K);
    Scal_grad_k=zeros(d,K);

    
    if iter==1 
        DkG_k=ones(1,K)*n/K*Bound_gamma;
    else
        DkG_k= ones(1,K)*norm(recv_grad*n/K);
    end
 
    
    
    %%% Optimal PA %%%%
    c_min=sqrt(P)*min(sqrt(pred_h')./kron(ones(1,T-iter+1),DkG_k'));  % K by T, max is conducted over differnet K.   

    if sum( (sqrt(2)*gamma)^2*c_min.^2/sigma_noise^2 )<res_dp
        c_t=c_min;
    else
        lb_sqrt_lambda=0;
        ub_sqrt_lambda=10^8;

        while ub_sqrt_lambda-lb_sqrt_lambda >10^(-8)

            mid_lambda_aircomp=(lb_sqrt_lambda+ub_sqrt_lambda)/2;

            if sum( (sqrt(2)*gamma)^2*  min( (1-mu/L).^-([iter:T]/2)/(sqrt(2*(1-mu/L)^(-T)*mid_lambda_aircomp)*gamma), c_min.^2/(sigma_noise)^2) ) <=res_dp
                ub_sqrt_lambda=mid_lambda_aircomp;
            else
                lb_sqrt_lambda=mid_lambda_aircomp;        
            end
        end
         mid_lambda_aircomp=(lb_sqrt_lambda+ub_sqrt_lambda)/2;
        c_t=min((1-mu/L)^(-iter/4)*sigma_noise*(2*(1-mu/L)^(-T)*mid_lambda_aircomp)^(-1/4)/sqrt(gamma), c_min);
    end  
    
   for k=1:K
        grad_k(:,k)=sum(grad_inv(sum(D(1:k))-D(k)+1:sum(D(1:k)),:),1)'/D(k)+ lambda*2*w_aircomp;
        
        Scal_grad_k(:,k)=grad_k(:,k)*D(k)*min(1, sqrt(P)/(c_t(1)/h_current(k))/norm(grad_k(:,k)*D(k)));
           
    end
    
    recv_grad=(sum(Scal_grad_k,2)+channel_noise(iter,:)'/c_t(1))/n;
    
    w_aircomp=w_aircomp-mu/L*recv_grad;
    
    %projection w
    if norm(w_aircomp)>max_w
        mark_aircomp=1;
        w_aircomp=w_aircomp/norm(w_aircomp)*max_w;
    end
    
%   record_w_aircomp=[record_w_aircomp,w_aircomp];
    
%     Loss_iter_aircomp(iter,NumUser_index,ave_index)=( sum((X*w_aircomp-Y).^2)/2/n+5*10^(-5)*norm(w_aircomp)^2 )/L_opt;
    
    
end  