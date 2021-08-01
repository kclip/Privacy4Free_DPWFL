clear all
clc
close all

load('Select Your Local Path/generated_Data_dim10.mat')

% data dimension
d = size(X,2); 

% num samples
n = size(X,1);



% estimate Optimal w
w_opt=inv(X'*X+5*10^(-5)*eye(d))*X'*Y;



lambda=5*10^(-4);





L_opt=sum((X*w_opt-Y).^2)/2/n+5*10^(-5)*norm(w_opt)^2;







%Channel noise, assume identical to all users
sigma_noise=1;  

%Power Constraint  
P=10^4;  

epsilon=20;

C=1.849;  % get from delta=0.01 
dp=(sqrt(epsilon+C^2)-C)^2;



Bound_gamma=20;  %clipped value
max_w=10;  % Limited searching space





% Channel Coefficients 
rho=0;
kappa=5;
 
 
K_set=[2 3 5 10 15 30]; % number of users 

T=30;  % Comm. Budget




Loss_aircomp=zeros(100,length(K_set));
Loss_BM_NA=zeros(100,length(K_set));
Loss_WithoutDP_NA=zeros(100,length(K_set));
Loss_TDMA=zeros(100,length(K_set));
Loss_BM_OA=zeros(100,length(K_set));
Loss_WithoutDP_OA=zeros(100,length(K_set));

rng(1)

for ave_index=1:100
    ave_index
    
    
    clear i
    
    gr_bais=sqrt(1/2)*(randn(T,max(K_set))+i*randn(T,max(K_set)));



    channel_noise=sigma_noise*randn(T,d);


   

    gr=zeros(T,max(K_set));

    gr(1,:)=gr_bais(1,:);

    for t=2:T
        gr(t,:)=rho*gr(t-1,:)+sqrt(1-rho^2)*gr_bais(t,:);
    end

    g=sqrt(kappa/(kappa+1))+ sqrt(1/(kappa+1))*gr;

    
    for K_index=1:6
        
        tic
        
        K=K_set(K_index); % comm. round

        % local data size
        D=[ones(1,K-1)*floor(n/K),n-floor(n/K)*(K-1)];
        
     
        h= abs(g(:,1:K));        
        
   
        TDMA_online
        Loss_TDMA(ave_index,K_index)= ( sum((X*w_TDMA-Y).^2)/2/n+5*10^(-5)*norm(w_TDMA)^2 -L_opt)/L_opt;

        Benchmark_Orthogonal_online
        Loss_BM_OA(ave_index,K_index)=( sum((X*w_OAbm-Y).^2)/2/n+5*10^(-5)*norm(w_OAbm)^2 -L_opt)/L_opt;

        WithoutDP_Orthogonal_online
        Loss_WithoutDP_OA(ave_index,K_index)=( sum((X*w_OA-Y).^2)/2/n+5*10^(-5)*norm(w_OA)^2 -L_opt)/L_opt;


        Aircomp_online
        Loss_aircomp(ave_index,K_index)=( sum((X*w_aircomp-Y).^2)/2/n+5*10^(-5)*norm(w_aircomp)^2 -L_opt)/L_opt;

        Benchmark_NonOrthogonal_online
         Loss_BM_NA(ave_index,K_index)=( sum((X*w_NAbm-Y).^2)/2/n+5*10^(-5)*norm(w_NAbm)^2 -L_opt )/L_opt;

        WithoutDP_NonOrthogonal_online
         Loss_WithoutDP_NA(ave_index,K_index)=( sum((X*w_NA-Y).^2)/2/n+5*10^(-5)*norm(w_NA)^2 -L_opt)/L_opt;


         toc

%         Aircomp_online_WellPredict







    end

semilogy(K_set,sum(Loss_aircomp,1)/ave_index,'LineWidth',2,'Color', [0    0.4470    0.7410])
hold on 
semilogy(K_set,sum(Loss_BM_NA,1)/ave_index,'--','LineWidth',2,'Color',[0    0.4470    0.7410])
hold on
semilogy(K_set,sum(Loss_WithoutDP_NA,1)/ave_index,':*','LineWidth',2,'Color',[0    0.4470    0.7410])
hold on


semilogy(K_set,sum(Loss_TDMA,1)/ave_index,'LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on
semilogy(K_set,sum(Loss_BM_OA,1)/ave_index,'--','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on
semilogy(K_set,sum(Loss_WithoutDP_OA,1)/ave_index,':*','LineWidth',2,'Color',[0.8500    0.3250    0.0980])

hold off

ylabel('Normalized Optimality Gap','fontsize',17)
xlabel('The Number of Users','fontsize',17)
h=legend('Adaptive PA','Static PA','Without DP constraint');


set(gca,'FontName','Times New Roman','FontSize',16);
legend boxoff 


set(h,'FontName','Times New Roman','FontSize',15)




pause(0.01)
end


