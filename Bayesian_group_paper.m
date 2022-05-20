function [Z,H]=Bayesian_group_paper(Y_all,X_all,N_point,etc,G)
Pm=[];
H=[];
rho=0.001;

[T,S,K]=size(Y_all);
N=size(X_all,2);
gap_d=1.085;
search_area=[-1:2/N_point:1]'/gap_d;
reslu=2/N_point/gap_d;
N_grid=length(search_area);
for kk=1:K
   search_area(:,kk)=[-1:2/N_point:1]';
   Fm(:,:,kk)=exp(-1i*pi*gap_d*(0:N-1)'*search_area(:,kk)');
end

a=0.0001;b=0.0001;
maxiter=500;
tol=1e-4;
%%%%%%%%%%%%%%%
delta_last=100;
converged = false;
iter = 0;

Aw=zeros(T,N_grid,K);
for kk=1:K
     Aw(:,:,kk)=X_all(:,:,kk)*Fm(:,:,kk);
end
Z=rand(K,G);
Z= diag(  1./sum(Z,2) ) *   Z;
Aw=[Aw,Aw];
%initialization
alpha0=1;
alpha_all=ones(N_grid*2, K);

alpha_all(N_grid+1:end,:)= alpha_all(N_grid+1:end,:)*rho;

mu=zeros(N_grid*2,S,K);
for kk=1:K
    Phi_delta = Aw(:,:,kk) *  diag(1./alpha_all(:,kk));
    V_temp= 1/alpha0*eye(T) + Phi_delta * Aw(:,:,kk)';
    Sigma(:,:,kk) = diag(1./alpha_all(:,kk)) -Phi_delta' * (V_temp \Phi_delta);
    mu(:,:,kk) = alpha0 * Sigma(:,:,kk) * Aw(:,:,kk)' * Y_all(:,:,kk);
end


alpha_star=rand(N_grid,G);
alpha_remain=ones(N_grid,K);

while ~converged

Phi=Aw;   
    
   switch iter-floor(iter/5)*5
       case 0               %update alpha0
               term0=zeros(T,S,K);
               term2=0;
               for kk=1:K
                   term0(:,:,kk)=Phi(:,:,kk)*mu(:,:,kk);
                   temp1=sum(diag(Phi(:,:,kk)* Sigma(:,:,kk)*Phi(:,:,kk)'));
                   term2=term2+temp1;                  
               end
               resid=Y_all-term0;
               alpha0=( S*K*T + a )/( b +  norm(resid(:), 'fro')^2  +   S* real(term2) );

       case 1               %update x
           
                  
               for kk=1:K
                   diag_inv= diag(  1./   (  [alpha_star  *  ( Z(kk,:).'  ) ; rho^(-1)* alpha_remain(:,kk)    ])  );
                   Phi_delta = Phi(:,:,kk) *  diag_inv;
                   V_temp= 1/alpha0*eye(T) + Phi_delta * Phi(:,:,kk)';
                   Sigma(:,:,kk) = diag_inv-Phi_delta' * (V_temp \Phi_delta);
                   mu(:,:,kk) =   alpha0 * Sigma(:,:,kk) * Phi(:,:,kk)' * Y_all(:,:,kk);
               end
               
               
       case 2               %update alpha^*   and  remained 
               mu2=abs(mu).^2;
               xx=mu2;
               for kk=1:K
                   sigma2_kk= real(   diag( Sigma(:,:,kk))  );
                   xx(:,:,kk)=mu2(:,:,kk) + sigma2_kk*ones(1,S);
               end
          
               c_k= a + sum(Z)*S;
               xx_temp=zeros(2*N_grid, K);
               for kk=1:K
                   xx_temp(:,kk)=  sum( xx(:,:,kk), 2);
               end
               d_ik= b +  xx_temp(1:N_grid,:)* Z;
   
               for gg=1:G
                    alpha_star(:,gg)=  c_k(gg)./d_ik(:,gg);
                    ln_alpha_star(:,gg)=   psi( c_k(gg) )   -log(d_ik(:,gg));
               end
               if iter<100
                 ln_alpha_star=log(alpha_star);
               end      
               
               
              %%%%%%%%%%%%%%% update alpha_remain
              alpha_remain =  (a+1)  ./   (  b+xx_temp(N_grid+1:end,:) *(rho^(-1))  );              

               

       case 3
             for kk=1:K
                 for gg=1:G
                    t1= S* sum( ln_alpha_star(:,gg)) ;
                    t2=  sum(xx_temp(1:N_grid,kk).* alpha_star(:,gg) )   ;            
                    et(kk,gg)= t1-t2;
                  end
             end
             Z1=[];
             
             for kk=1:K
                 temp(kk,:)= exp(et(kk,:)-max(max(et(kk,:)))); 
             end
  
             Z= diag(  1./sum(temp,2) ) *   temp;


             
       case 4                
              %%%%%%%%%%%%%%%%%%%%%%%%%%%% grid refine
              for kk=1:K
                    term0(:,:,kk)=Phi(:,:,kk)*mu(:,:,kk);
              end
              resid=Y_all-term0;   
              m1=mu(1:N_grid,:,:); m2=mu(N_grid+1:end,:,:);
              mu12=m1+m2;
              mu222=  abs(mu12).^2;
            
              for kk=1:K

                        sum_mu=sum(  mu222(:,:,kk) , 2);
                        Pm=sum_mu;
                        [~,sort_ind]=sort(Pm, 'descend');    
                        index_amp = sort_ind(1:etc);
                        
                        Sigma12= Sigma(1:N_grid,1:N_grid,kk)+Sigma(N_grid+1:end,N_grid+1:end,kk)+ Sigma(1:N_grid,N_grid+1:end,kk) +  Sigma(N_grid+1:end,1:N_grid,kk);
                        
                        tempPS=Phi(:,1:N_grid,kk)*  Sigma12(:,index_amp) ;    
                        df=zeros(length(index_amp),1);
                        for j=1:length(index_amp)
                                ii=index_amp(j);
                                ai=Fm(:,ii,kk);
                                mut=mu12(ii,:,kk);
                                Sigmat=Sigma12(:,ii);
                                c1=mut*mut' +  S*Sigmat(ii);
                                c1=abs(c1)*(-alpha0);
                                Yti=resid(:,:,kk) +  Phi(:, ii,kk)*mu12(ii,:,kk);
                                c2= S*(  tempPS(:,j) - Phi(:,ii,kk)*Sigmat(ii) )  -Yti*(mut');
                                c2= c2*(-alpha0);
                                phii=Phi(:,ii,kk);
                                sinta= search_area(ii,kk); costa=cos(asin(sinta));
                                c3=(-1i*pi*gap_d* costa)*[0:N-1]';    %  c3=(-1i*2*pi*gap_d/sqrt(N))*[0:N-1]';
                                tt1=  X_all(:,:,kk)*(c3.*ai); 
                                f1= tt1'*phii*c1  +   tt1'*c2;
                                f1= 2*real(f1);
                                df(j)=f1;
%                                 angle_cand= sign(f1)*reslu/100 ;       
%                                 sin_add = search_area(ii,kk) + angle_cand;
%                                 search_area(ii,kk)=sin_add;
%                                 ro1=exp(-sin_add*pi*gap_d*1i);            
%                                 Fm(:,ii,kk)=ro1.^((0:N-1)')/sqrt(N);
%                                 Aw(:,ii,kk)=X_all(:,:,kk)* Fm(:,ii,kk);
%                                 Aw(:,N_grid+ii,kk)= Aw(:,ii,kk);
                        end   
                        ddff=sign(df)*reslu/100;    
                        search_area(index_amp,kk) = search_area(index_amp,kk) +  ddff;
                        Fm(:,index_amp,kk)=exp(-1i*pi*gap_d*(0:N-1)'* (search_area(index_amp,kk))');
                        Aw(:,index_amp,kk)=X_all(:,:,kk)* Fm(:,index_amp,kk);
                        Aw(:,N_grid+index_amp,kk)=Aw(:,index_amp,kk);
              end
 
   end
   

    % stopping criteria
      erro=norm(alpha_all - delta_last)/norm(delta_last);
    if erro < tol || iter >= maxiter
        converged = true;
    end
     iter = iter + 1;
   
end

% figure;hold on;
% for kk=1:K
%    subplot(4,5,kk);plot(sum(mu2(:,:,kk),2));   title(kk);    
% end


for kk=1:K
    Pm=sum(mu222(:,:,kk),2);
    [~,sort_ind]=sort(Pm, 'descend');  
    ther1=mean(Pm(sort_ind(1:etc)))*0.01;   
    index_amp=find(Pm>ther1);
    Tn= X_all(:,:,kk)*Fm(:,index_amp,kk);
    Ss= Tn \ Y_all(:,:,kk);
    H(:,:,kk)= Fm(:,index_amp,kk)*Ss;    
end

%  figure(5);plot(1./alpha_star)
%  figure(6);plot(1./alpha_remain)         
%                 
               










