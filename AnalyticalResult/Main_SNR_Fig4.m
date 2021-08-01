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

% local data size
D=n/K;

%regularization weight
lambda=5*10^(-5);

eigval=svd(X'*X/n+2*lambda*eye(d));
L=max(eigval);
mu=min(eigval);

% estimate Optimal w
w_opt=inv(X'*X+2*n*lambda*eye(d))*X'*Y;

% Limited searching space
max_w=norm(w_opt);

% %estimate bounded sample-wise gradient from the collected data 
% beta2=max(sum( ((Y-X*w_opt).*(-X)).^2, 2 ));
% ave_g_norm0= norm(sum(Y.*(-X),1)/n);
% max_g_norm0=  max( sqrt( sum( ( Y.*(-X) ).^2 ,2) ) );
% beta1=(max_g_norm0^2-beta2)/ave_g_norm0^2;


%estimate bound of gradient 
for i=1:K
    L_invK=0;
    L_max_invK=[];
    for j=1:D
        L_invK=L_invK+ X(D*(i-1)+j,:)'*X(D*(i-1)+j,:);
        L_max_invK=[L_max_invK, max(svd( X(D*(i-1)+j,:)'*X(D*(i-1)+j,:) + 2*lambda*eye(d) ))];
    end        
   Bound_G(i)=max(svd(L_invK/D+2*lambda*eye(d)))*2*max_w;
   Bound_gamma(i)=max(L_max_invK)*2*max_w;
end
Bound_gamma=max(Bound_gamma);




L_opt=sum((X*w_opt-Y).^2)/2/n+lambda*norm(w_opt)^2;

L_1=sum((-Y).^2)/2/n;

%Initial distance F(w1)-F(w*)=F0
F0=L_1-L_opt;


epsilon=20;
delta=0.01;
C= ((1/delta)*exp(-lambertw(0, (2*(1/delta)^2)/pi)/2))/sqrt(pi);
% dp_analog=(sqrt(epsilon_tot+C_analog^2)-C_analog)^2;
% C=1.849;  % get from delta=0.01 
dp=(sqrt(epsilon+C^2)-C)^2;

D=n/K*ones(1,K);


T=30;


%Channel noise, assume identical to all users
sigma_noise=1;  

OptGap_aircomp=zeros(1000,41);
OptGap_BM_NA=zeros(1000,41);
OptGap_WithoutDP_NA=zeros(1000,41);
OptGap_TDMA=zeros(1000,41);
OptGap_BM_OA=zeros(1000,41);
OptGap_WithoutDP_OA=zeros(1000,41);



for ave_index=1:1000
    
    kappa=10;
    clear i
    g=sqrt(kappa/(kappa+1))+sqrt(1/(kappa+1))*sqrt(1/2)*(randn(1,K)+i*randn(1,K));
    h=kron(ones(T,1),abs(g));

%     h=ones(T,K);

    for snr_index=1:41



    %     rng(1)
        %Channel realization for single user  
    %     h=exprnd(ones(T,K)); 
         
    %      h=abs(randn(T,K));
    %     %large_scale=[1:dis_index:dis_index*K]/K;
    %     large_scale=[1,(1-mu/L)^(T/2)];
    %     h=kron(ones(T,1),large_scale);
        %h=h_0.*kron(ones(T,1),large_scale);

        %Power Constraint 
        P=10^(((snr_index-1+5)/10))*sigma_noise^2*d;

     



        TDMA_analytical
        OptGap_TDMA(ave_index,snr_index)=(1-mu/L)^(T/K)*( F0+d/2/L/n^2* ( (1-mu/L).^(-[1:T/K])*sum(x_TDMA,2)) )/L_opt;
        OptGap_WithoutDP_OA(ave_index,snr_index)=(1-mu/L)^(T/K)*( F0+d/2/L/n^2* ( (1-mu/L).^(-[1:T/K])*sum(x_WithoutDP_OA,2)) )/L_opt;

        Aircomp_analytical
        OptGap_aircomp(ave_index,snr_index)=(1-mu/L)^T*( F0 + d/2/L/n^2* ( (1-mu/L).^(-[1:T])*x_aircomp') )/L_opt;
        OptGap_WithoutDP_NA(ave_index,snr_index)=(1-mu/L)^T*( F0 + d/2/L/n^2* ( (1-mu/L).^(-[1:T])*x_WithoutDP_NA') )/L_opt;

        Benchmark_Orthogonal
        OptGap_BM_OA(ave_index,snr_index)=(1-mu/L)^(T/K)*( F0+d/2/L/n^2* ( (1-mu/L).^(-[1:T/K])*sum(x_BM_oa,2)) )/L_opt;

        Benchmark_NonOrthogonal
        OptGap_BM_NA(ave_index,snr_index)=(1-mu/L)^T*( F0+d/2/L/n^2* ( (1-mu/L).^(-[1:T])*x_BM_na') )/L_opt;

    end

% semilogy( ([1:11]-1)*2,OptGap_aircomp,'LineWidth',2)
% hold on 
% semilogy( ([1:11]-1)*2,OptGap_TDMA,'--','LineWidth',2)
% hold on 
% semilogy( ([1:11]-1)*2,OptGap_BM_OA,'-.','LineWidth',2)
% hold on 
% semilogy( ([1:11]-1)*2,OptGap_BM_NA,':','LineWidth',2)
snr_index=[1:41];
P=10.^(((snr_index-1+5)/10))*sigma_noise^2*d;

snr_range=10*log10(P/(d*sigma_noise^2));


% for plot the legend
semilogy(snr_range,sum(OptGap_aircomp,1)/ave_index,'LineWidth',2,'Color', 'black')
hold on 
semilogy(snr_range,sum(OptGap_BM_NA,1)/ave_index,'--','LineWidth',2,'Color','black')
hold on 
semilogy(snr_range,sum(OptGap_WithoutDP_NA,1)/ave_index,':','LineWidth',2,'Color','black')
hold on 

% for plot NOMA
semilogy(snr_range,sum(OptGap_aircomp,1)/ave_index,'LineWidth',2,'Color', [0    0.4470    0.7410])
hold on 
semilogy(snr_range,sum(OptGap_BM_NA,1)/ave_index,'--','LineWidth',2,'Color',[0    0.4470    0.7410])
hold on 
semilogy(snr_range,sum(OptGap_WithoutDP_NA,1)/ave_index,':','LineWidth',2,'Color',[0    0.4470    0.7410])
hold on 


% for plot OMA
semilogy(snr_range,sum(OptGap_TDMA,1)/ave_index,'LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on 
semilogy(snr_range,sum(OptGap_BM_OA,1)/ave_index,'--','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on 
semilogy(snr_range,sum(OptGap_WithoutDP_OA,1)/ave_index,':','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold off

pause(0.1)

end

figure(1)
hold on

index_NA=find((sum(OptGap_aircomp,1)/ave_index-sum(OptGap_WithoutDP_NA,1)/ave_index)<10^(-3));
scatter(snr_range(index_NA(length(index_NA))),sum(OptGap_aircomp(:,index_NA(length(index_NA))))/ave_index,50,'filled','MarkerFaceColor',[0    0.4470    0.7410])

index_OA=find((sum(OptGap_TDMA,1)/ave_index-sum(OptGap_WithoutDP_OA,1)/ave_index)<10^(-3));
scatter(snr_range(index_OA(length(index_OA))),sum(OptGap_TDMA(:,index_OA(length(index_OA))))/ave_index,50,'filled','MarkerFaceColor',[0.8500    0.3250    0.0980])


% plot([1:10]*10,(OptGap_aircomp),'LineWidth',2)
% hold on 
% plot([1:10]*10,(OptGap_TDMA),'--','LineWidth',2)
% hold on 
% plot([1:10]*10,(OptGap_BM_OA),'-.','LineWidth',2)
% hold on 
% plot([1:10]*10,(OptGap_BM_NA),':','LineWidth',2)
ylabel('Normalized Optimality Gap','fontsize',17)
xlabel('SNR (dB)','fontsize',17)
%  xticks([-10:2:40])
% legend('Non-orthogonal Aggregation','Orthogonal Access','Benchmark Orthogonal','Benchmark Non-orthogonal')
% legend('Adaptive PA','Constant PA','Without DP constraint');

 h=legend('Adaptive PA','Static PA','Without DP constraint');


    set(gca,'FontName','Times New Roman','FontSize',16);
legend boxoff 


set(h,'FontName','Times New Roman','FontSize',15)
ylim([5*10^-3,1.5*10])
xlim([12,38])