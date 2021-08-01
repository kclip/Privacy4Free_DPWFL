clear all
clc
close all

load('Select Your Local Path/generated_Data_dim10.mat')

% data dimension
d = size(X,2); 

% num samples
n = size(X,1);

% num user
K=10;

%regularization weight
lambda=5*10^(-5);

% estimate Optimal w
% wrong w_opt=inv(X'*X+5*10^(-5)*eye(d))*X'*Y;
w_opt=inv(X'*X+2*n*lambda*eye(d))*X'*Y;

% Limited searching space
max_w=norm(w_opt);

L_opt=sum((X*w_opt-Y).^2)/2/n+lambda*norm(w_opt)^2;

L_1=sum((-Y).^2)/2/n;

%Initial distance F(w1)-F(w*)=F0
F0=L_1-L_opt;


%comm. round
T=30;
    
%DP setting
epsilon=20;
C=1.849;  % get from delta=0.01 
dp=(sqrt(epsilon+C^2)-C)^2;


%Channel Noise 
sigma_noise=1; 

%Power Constraint 
P=10^4;


norm_X=[[1:n]',sum(X.^2,2)];
rank_index=sortrows(norm_X,2);



Bound_gamma=max(sum(X.^2,2)+2*lambda)*2*max_w;

OptGap_aircomp=zeros(1000,20);
OptGap_BM_NA=zeros(1000,20);
OptGap_WithoutDP_NA=zeros(1000,20);
OptGap_TDMA=zeros(1000,20);
OptGap_BM_OA=zeros(1000,20);
OptGap_WithoutDP_OA=zeros(1000,20);
          

rng(1)

for ave_index=1:1000
    
    ave_index
    kappa=10;
    clear i
    g=sqrt(kappa/(kappa+1))+sqrt(1/(kappa+1))*sqrt(1/2)*(randn(1,K)+i*randn(1,K));
    h=kron(ones(T,1),abs(g));


    for het_index=1:20



            max_ratio=1/K+(het_index-1)*(1-1/K)/20;

            other_indvNum=round((1-max_ratio)/(K-1)*n);

            % num D_k for all devices K, vector 1 by K
            D=[n-other_indvNum*(K-1), other_indvNum*ones(1,K-1)];

            distr_index=rank_index([n/2-other_indvNum/2*K+1:n/2+other_indvNum/2*K],1);
            distr_index=distr_index(randperm(length(distr_index)));


            %For the different distribution of D_k, G_k may vary a lot 
            %The above process is trying to adjust G_k almost same for all k

            % index for calculation of G_k 
            calG_index=[rank_index([[1:n/2-other_indvNum/2*K],[n/2+other_indvNum/2*K+1:n]],1);distr_index];

            for i=1:K
                L_invK=0;           
                for j=0:D(i)-1   
                    L_invK=L_invK+ X(calG_index(sum(D(1:i))-j),:)'* X(calG_index(sum(D(1:i))-j),:);  
    %                 L_invK=L_invK+ X(sum(D(1:i))-j,:)'* X(sum(D(1:i))-j,:);   
                end        
                Bound_G(i)=max(svd(L_invK/D(i)+10^(-4)*eye(d)))*2*max_w;

            end

            TDMA_analytical
            OptGap_TDMA(ave_index,het_index)=(1-mu/L)^(T/K)*( F0+d/2/L/n^2* ( (1-mu/L).^(-[1:T/K])*sum(x_TDMA,2)) )/L_opt;
            OptGap_WithoutDP_OA(ave_index,het_index)=(1-mu/L)^(T/K)*( F0+d/2/L/n^2* ( (1-mu/L).^(-[1:T/K])*sum(x_WithoutDP_OA,2)) )/L_opt;

            Aircomp_analytical
            OptGap_aircomp(ave_index,het_index)=(1-mu/L)^T*( F0 + d/2/L/n^2* ( (1-mu/L).^(-[1:T])*x_aircomp') )/L_opt;
            OptGap_WithoutDP_NA(ave_index,het_index)=(1-mu/L)^T*( F0 + d/2/L/n^2* ( (1-mu/L).^(-[1:T])*x_WithoutDP_NA') )/L_opt;


            Benchmark_Orthogonal
            OptGap_BM_OA(ave_index,het_index)=(1-mu/L)^(T/K)*( F0+d/2/L/n^2* ( (1-mu/L).^(-[1:T/K])*sum(x_BM_oa,2)) )/L_opt;

            Benchmark_NonOrthogonal
            OptGap_BM_NA(ave_index,het_index)=(1-mu/L)^T*( F0+d/2/L/n^2* ( (1-mu/L).^(-[1:T])*x_BM_na') )/L_opt;

    end


semilogy(1/K+([1:20]-1)*(1-1/K)/20,sum(OptGap_aircomp,1)/ave_index,'LineWidth',2,'Color', [0    0.4470    0.7410])
hold on 
semilogy(1/K+([1:20]-1)*(1-1/K)/20,sum(OptGap_BM_NA,1)/ave_index,'--','LineWidth',2,'Color',[0    0.4470    0.7410])
hold on 
semilogy(1/K+([1:20]-1)*(1-1/K)/20,sum(OptGap_WithoutDP_NA,1)/ave_index,':','LineWidth',2,'Color',[0    0.4470    0.7410])
hold on 

semilogy(1/K+([1:20]-1)*(1-1/K)/20,sum(OptGap_TDMA,1)/ave_index,'LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on 
semilogy(1/K+([1:20]-1)*(1-1/K)/20,sum(OptGap_BM_OA,1)/ave_index,'--','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on 
semilogy(1/K+([1:20]-1)*(1-1/K)/20,sum(OptGap_WithoutDP_OA,1)/ave_index,':','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold off 

pause(0.01)
end
figure(1)
hold on 
ylabel('Normalized Optimality Gap','fontsize',17)
xlabel('max D_k/D_{tot}','fontsize',17)
h=legend('Adaptive PA','Static PA','Without DP constraint');

xlim([0.1,0.95])

set(gca,'FontName','Times New Roman','FontSize',16);
legend boxoff 


set(h,'FontName','Times New Roman','FontSize',15)