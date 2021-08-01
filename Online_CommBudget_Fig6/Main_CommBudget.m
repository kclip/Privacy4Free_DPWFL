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

L_1=sum((-Y).^2)/2/n;

%Initial distance F(w1)-F(w*)=F0
F0=L_1-L_opt;




%Channel noise, assume identical to all users
sigma_noise=1;  

%Power Constraint  
P=10^4;  

epsilon=20;

C=1.849;  % get from delta=0.01 
dp=(sqrt(epsilon+C^2)-C)^2;



Bound_gamma=20;  %clipped value
max_w=10;  % Limited searching space





rng(1)
rho=0;
K=10; % number of users

T_set=[1 2 3 5 10 20 30 40 50 60 70];   %set of comm.budget 

for ave_index=1:100
    ave_index
    
    
    clear i
    
    gr_bais=sqrt(1/2)*(randn(max(T_set),K)+i*randn(max(T_set),K));



    channel_noise_set=sigma_noise*randn(max(T_set),d);


    kappa=5;

    gr=zeros(max(T_set),K);

    gr(1,:)=gr_bais(1,:);

    for t=2:max(T_set)
        gr(t,:)=rho*gr(t-1,:)+sqrt(1-rho^2)*gr_bais(t,:);
    end

    g=sqrt(kappa/(kappa+1))+ sqrt(1/(kappa+1))*gr;

    
    for T_index=1:11
        
        tic
        
        T=T_set(T_index); % comm. round

        % local data size
        D=[ones(1,K-1)*floor(n/K),n-floor(n/K)*(K-1)];
        
     
        h= abs(g(1:T,:));
        channel_noise=channel_noise_set(1:T,:);
        
        
        Loss_initial=( sum((-Y).^2)/2/n )/L_opt;
        if T<10
            Aircomp_online
            Loss_aircomp(ave_index,T_index)=( sum((X*w_aircomp-Y).^2)/2/n+5*10^(-5)*norm(w_aircomp)^2 -L_opt)/L_opt;

            Benchmark_NonOrthogonal_online
            Loss_BM_NA(ave_index,T_index)=( sum((X*w_NAbm-Y).^2)/2/n+5*10^(-5)*norm(w_NAbm)^2 -L_opt )/L_opt;

            WithoutDP_NonOrthogonal_online
            Loss_WithoutDP_NA(ave_index,T_index)=( sum((X*w_NA-Y).^2)/2/n+5*10^(-5)*norm(w_NA)^2 -L_opt)/L_opt;
        else
            TDMA_online
            Loss_TDMA(ave_index,T_index-4)= ( sum((X*w_TDMA-Y).^2)/2/n+5*10^(-5)*norm(w_TDMA)^2 -L_opt)/L_opt;

            Benchmark_Orthogonal_online
            Loss_BM_OA(ave_index,T_index-4)=( sum((X*w_OAbm-Y).^2)/2/n+5*10^(-5)*norm(w_OAbm)^2 -L_opt)/L_opt;

            WithoutDP_Orthogonal_online
            Loss_WithoutDP_OA(ave_index,T_index-4)=( sum((X*w_OA-Y).^2)/2/n+5*10^(-5)*norm(w_OA)^2 -L_opt)/L_opt;


            Aircomp_online
            Loss_aircomp(ave_index,T_index)=( sum((X*w_aircomp-Y).^2)/2/n+5*10^(-5)*norm(w_aircomp)^2 -L_opt)/L_opt;

            Benchmark_NonOrthogonal_online
             Loss_BM_NA(ave_index,T_index)=( sum((X*w_NAbm-Y).^2)/2/n+5*10^(-5)*norm(w_NAbm)^2 -L_opt )/L_opt;

            WithoutDP_NonOrthogonal_online
             Loss_WithoutDP_NA(ave_index,T_index)=( sum((X*w_NA-Y).^2)/2/n+5*10^(-5)*norm(w_NA)^2 -L_opt)/L_opt;
        end

         toc

%         Aircomp_online_WellPredict







    end

semilogy(T_set,sum(Loss_aircomp,1)/ave_index,'LineWidth',2,'Color', [0    0.4470    0.7410])
hold on 
semilogy(T_set,sum(Loss_BM_NA,1)/ave_index,'--','LineWidth',2,'Color',[0    0.4470    0.7410])
hold on
semilogy(T_set,sum(Loss_WithoutDP_NA,1)/ave_index,':*','LineWidth',2,'Color',[0    0.4470    0.7410])
hold on


semilogy([10 20 30 40 50 60 70],sum(Loss_TDMA,1)/ave_index,'LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on
semilogy([10 20 30 40 50 60 70],sum(Loss_BM_OA,1)/ave_index,'--','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on
semilogy([10 20 30 40 50 60 70],sum(Loss_WithoutDP_OA,1)/ave_index,':*','LineWidth',2,'Color',[0.8500    0.3250    0.0980])

hold off

ylabel('Normalized Optimality Gap','fontsize',17)
xlabel('Communication Budget','fontsize',17)
h=legend('Adaptive PA','Static PA','Without DP constraint');


set(gca,'FontName','Times New Roman','FontSize',16);
legend boxoff 


set(h,'FontName','Times New Roman','FontSize',15)





pause(0.01)
end


