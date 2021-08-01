x_BM_oa=zeros(T/K,K);

for k=1:K


min_x_oa=(D(k)*Bound_G(k))^2*(sigma_noise./h([k:K:T],k)).^2/P;

x_BM_oa(:,k)=max(  (sqrt(2)*Bound_gamma)^2/ ( dp/(T/k) ) ,  min_x_oa);


end