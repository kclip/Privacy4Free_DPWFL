clear all
clc
close all

load('Select Your Local Path/MNIST.mat')



X= train_set(:, 2:size(train_set,2));
% num samples
n = size(X,1);

Y= zeros(n,10);
i = [1:n]'; 
j = train_set(:,1)+1;
Y(sub2ind(size(Y), i, j))=1;
X=[ones(n,1),X];
clear j

X_test=[ ones(size(test_set,1),1), test_set(:,2:785)]; 
Y_test=test_set(:,1);

% data dimension
d = size(X,2); 

%model dimension, include bias  
dm=d*10; 

%Weight for regularization 
lambda=10^-2;


% LogLoss:  -sum(sum(Y.*log( exp(X*W_aircomp)./sum(exp(X*W_aircomp),2) )))/n+lambda*sum(sum(W_aircomp.^2));

% H_matrix=kron(1/10*eye(10)-1/10/10*ones(10,10),X'*X)/n + 2*lambda*eye(dm);
% L=max(svd(H_matrix));
% mu=2*lambda;
% 
% clear H_matrix 

mu=2*lambda; 
L=2.5;
clear train_set test_set


%Channel noise, assume identical to all users
sigma_noise=1;  

%Power Constraint  
P=10^4;  

epsilon=2;
delta_tot=1e-2;
C= ((1/delta_tot)*exp(-lambertw(0, (2*(1/delta_tot)^2)/pi)/2))/sqrt(pi);
dp=(sqrt(epsilon+C^2)-C)^2;



Bound_gamma=50;  %clipped value
max_w=10;  % Limited searching space




%learning rate
eta=1/L; 

rho=0;
K=10; % number of users

% local data size
D=[ones(1,K-1)*floor(n/K),n-floor(n/K)*(K-1)];

T_set=[10 30 60 100 150 220 300];   %set of comm.budget 

for ave_index=1:100
    ave_index
    
    W0=zeros(d, 10);
    
%     rnd_index = randperm(n);
%     X=X(rnd_index,:);
%     Y=Y(rnd_index,:);
   
    clear i
    
    gr_bais=sqrt(1/2)*(randn(max(T_set),K)+i*randn(max(T_set),K));



    channel_noise_set=sigma_noise*randn(max(T_set),dm);


    kappa=5;

    gr=zeros(max(T_set),K);

    gr(1,:)=gr_bais(1,:);

    for t=2:max(T_set)
        gr(t,:)=rho*gr(t-1,:)+sqrt(1-rho^2)*gr_bais(t,:);
    end

    g=sqrt(kappa/(kappa+1))+ sqrt(1/(kappa+1))*gr;

    
    for T_index=1:length(T_set)
        
        tic
        
        T=T_set(T_index); % comm. round

       
     
        h= abs(g(1:T,:));
        channel_noise=channel_noise_set(1:T,:);
        
        
        

            TDMA_online_Softmax

            Loss_TDMA(ave_index,T_index)=-sum(sum(Y.*log( exp(X*W_TDMA)./sum(exp(X*W_TDMA),2) )))/n+lambda*sum(sum(W_TDMA.^2));
            [M,pred_label]= max(exp(X_test*W_TDMA)./sum(exp(X_test*W_TDMA),2),[],2);
            pred_label=pred_label-1;    
            Err_TDMA(ave_index,T_index)= sum ( Y_test~= pred_label) /length(Y_test);



            Benchmark_Orthogonal_online_Softmax

            Loss_BM_OA(ave_index,T_index)=-sum(sum(Y.*log( exp(X*W_OAbm)./sum(exp(X*W_OAbm),2) )))/n+lambda*sum(sum(W_OAbm.^2));
            [M,pred_label]= max(exp(X_test*W_OAbm)./sum(exp(X_test*W_OAbm),2),[],2);
            pred_label=pred_label-1;    
            Err_BM_OA(ave_index,T_index)= sum ( Y_test~= pred_label) /length(Y_test);



            WithoutDP_Orthogonal_online_Softmax
            Loss_WithoutDP_OA(ave_index,T_index)=-sum(sum(Y.*log( exp(X*W_OA)./sum(exp(X*W_OA),2) )))/n+lambda*sum(sum(W_OA.^2));                       
            [M,pred_label]= max(exp(X_test*W_OA)./sum(exp(X_test*W_OA),2),[],2);
            pred_label=pred_label-1;    
            Err_WithoutDP_OA(ave_index,T_index)= sum ( Y_test~= pred_label) /length(Y_test);


        

         toc
         
         end

   






figure(1)

plot(T_set,sum(Loss_TDMA,1)/ave_index,'LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on
plot(T_set,sum(Loss_BM_OA,1)/ave_index,'--','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on
plot(T_set,sum(Loss_WithoutDP_OA,1)/ave_index,':*','LineWidth',2,'Color',[0.8500    0.3250    0.0980])

hold off

ylabel('Traning Loss','fontsize',17)
xlabel('Communication Budget','fontsize',17)
h=legend('Adaptive PA','Static PA','Without DP constraint');


set(gca,'FontName','Times New Roman','FontSize',16);
legend boxoff 


set(h,'FontName','Times New Roman','FontSize',15)

pause(0.01)


figure(2)


plot(T_set,sum(Err_TDMA,1)/ave_index,'LineWidth',2,'Color', [0.8500    0.3250    0.0980])
hold on
plot(T_set,sum(Err_BM_OA,1)/ave_index,'--','LineWidth',2,'Color',[0.8500    0.3250    0.0980])
hold on
plot(T_set,sum(Err_WithoutDP_OA,1)/ave_index,':*','LineWidth',2,'Color',[0.8500    0.3250    0.0980])


hold off

ylabel('Test Error','fontsize',17)
xlabel('Communication Budget','fontsize',17)
h=legend('Adaptive PA','Static PA','Without DP constraint');


set(gca,'FontName','Times New Roman','FontSize',16);
legend boxoff 


set(h,'FontName','Times New Roman','FontSize',15)

pause(0.01)

end


