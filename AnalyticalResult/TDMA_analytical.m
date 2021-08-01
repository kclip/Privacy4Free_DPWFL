x_TDMA=zeros(T/K,K);
x_WithoutDP_OA=zeros(T/K,K);

for k=1:K

min_x=(D(k)*Bound_G(k))^2*(sigma_noise./h([k:K:T],k)).^2/P;

x_WithoutDP_OA(:,k)=min_x;


if sum( (sqrt(2)*Bound_gamma)^2./min_x)<dp
    sqrt_lambda=0;
    x_TDMA(:,k)=min_x;
else
    lb_sqrt_lambda=0;
    ub_sqrt_lambda=10^9;
%     sum( Bound_gamma(k)^2./ ( max(ub_sqrt_lambda*(1-mu/L).^([1:T/K]'/2)*Bound_gamma(k), min_x) ))-dp;
    
    while ub_sqrt_lambda-lb_sqrt_lambda >10^(-9)
        
        mid_lambda=(lb_sqrt_lambda+ub_sqrt_lambda)/2;
        
        if sum( (sqrt(2)*Bound_gamma).^2./ ( max(sqrt(2)*mid_lambda*(1-mu/L)^(-T/K/2)*(1-mu/L).^([1:T/K]'/2)*Bound_gamma, min_x) ))<= dp
            ub_sqrt_lambda=mid_lambda;
        else
            lb_sqrt_lambda=mid_lambda;        
        end
    end
     mid_lambda=(lb_sqrt_lambda+ub_sqrt_lambda)/2;
    x_TDMA(:,k)=max(sqrt(2)*mid_lambda*(1-mu/L)^(-T/K/2)*(1-mu/L).^([1:T/K]'/2)*Bound_gamma, min_x);
end


end

% mid_lambda
% dp