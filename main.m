clear;

K=30;        % number of users
G_true=3;         % number of groups
N=80;        % number of antennas
T=60;        % Snapshot
N_point=N;   % number of grid points
SNR=10;
sigm=5;
L_common=2;
L_indiv=1;
L=L_common+L_indiv;
%% Generate the received signal

Ps=sqrt((   10.^(SNR/10)   )/2);
X=Ps*( randn(T,N)+1i*randn(T,N)  );

DoAs=zeros(G_true,L_common);
etc=10*(L);  
for gg=1:G_true
      DoAs(gg,:)= (rand(1,L_common)-0.5)*2*90;   
end
DOA_all=zeros(K,L);
for kk=1:K
    g_ind=randperm(G_true,1);
    randn_doa= (rand(1,L_indiv)-0.5)*2*90; 
    DOA_all(kk,:)=[DoAs(g_ind,:), randn_doa];
    Z(kk,1)=g_ind;
 end 
 for kk=1:K   
    aod_intro=DOA_all(kk,:);
    aodss=ones(10,1)*aod_intro + (rand(10,length(aod_intro))-0.5)*2*sigm;
    aaaa=zeros(1,length(aod_intro),10); 
    for iii=1:10      
       aaaa(1,:,iii)= aodss(iii,:);
    end 
    DOA_kk=aaaa(:)';
    A = exp(-1i*pi*(0:N-1)'*sind(DOA_kk));
    H=A*(randn(L*10,1)+1i*randn(L*10,1))/10;  
    noise=sqrt(1/2)*(randn(T,1)+1i*randn(T,1));
    H_all(:,:,kk)=H;    
    Y_all(:,:,kk)=X*H + noise;
    X_all(:,:,kk)=X;
 end  


%% The proposd method
[Z_our,H_est_our]=Bayesian_group_paper(Y_all,X_all,N_point,etc,G_true);
norm(H_all(:)-H_est_our(:),'fro')^2/norm(H_all(:))^2

