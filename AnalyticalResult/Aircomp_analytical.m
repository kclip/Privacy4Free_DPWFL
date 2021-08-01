% x_aircomp=zeros(T,K);


min_x=sigma_noise^2/P* max( ( kron(ones(T,1),D.^2.*Bound_G.^2)./(h.^2) )' ) ;  % K by T, max is conducted over differnet K.   
x_WithoutDP_NA=min_x;


if sum( (sqrt(2)*Bound_gamma)^2./min_x)<dp
    sqrt_lambda=0;
    x_aircomp=min_x;
else
    lb_sqrt_lambda=0;
    ub_sqrt_lambda=10^9;
%     sum( Bound_gamma(k)^2./ ( max(ub_sqrt_lambda*(1-mu/L).^([1:T/K]'/2)*Bound_gamma(k), min_x) ))-dp;
    
    while ub_sqrt_lambda-lb_sqrt_lambda >10^(-9)
        
        mid_lambda_aircomp=(lb_sqrt_lambda+ub_sqrt_lambda)/2;
        
        if sum( (sqrt(2)*Bound_gamma)^2./ ( max( sqrt(2)*mid_lambda_aircomp*(1-mu/L)^(-T/2)*(1-mu/L).^([1:T]/2)*Bound_gamma, min_x ) ))<= dp
            ub_sqrt_lambda=mid_lambda_aircomp;
        else
            lb_sqrt_lambda=mid_lambda_aircomp;        
        end    
   
    end
    
    mid_lambda_aircomp=(lb_sqrt_lambda+ub_sqrt_lambda)/2;
    x_aircomp=  max( sqrt(2)*mid_lambda_aircomp*(1-mu/L)^(-T/2)*(1-mu/L).^([1:T]/2)*Bound_gamma, min_x );
    
end

%      diff_aircomp=x_aircomp-min_x

% mid_lambda_aircomp
% 
% dp