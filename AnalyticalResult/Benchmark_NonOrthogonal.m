x_BM_na=zeros(T,K);


min_x_na=sigma_noise^2/P* max( ( kron(ones(T,1),D.^2.*Bound_G.^2)./(h.^2) )' );


x_BM_na =max(  (sqrt(2)*Bound_gamma)^2/ ( dp/T ) ,  min_x_na);


