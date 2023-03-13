%%
%%    By Hirokuni Iiboshi (2016)
%% Å@"Monetary Policy Regime Shifts under the Zero Lower Bound:
%%    An application of a stochastic rational expectations equilibrium to a Markov switching DSGE model" 
%%    Economic Modelling, Vol. 52, p186-205  
%%
%%

addpath('functions/');


%% save files of policy functions  
flag_save_policy_func = 1;

%%   Option of ZLB
const_ZLB = 1;   %  1-> ZLB,  0-> non ZLB

forward_looking = 1 % 1: perfect foresight, 2: stochastic  

% setting ranges of structural shocks
mlow_u = -0.1; mup_u = 0.1;  % markup shock
N_u = round((mup_u - mlow_u)/0.005);

mlow_g = -0.4; mup_g = 0.4; % real rate shock
N_g = round((mup_g - mlow_g)/0.002);


%% setting parameters based on Adam and Billi (2007, JME)
beta = 0.9913 ;  % (1+0.035/4)^(-1)
r_star = 1/beta -1 ;
sigma = 6.25;
kappa = 0.024;
lambda = 0.048/4^2 ; % 0.003
rho_u = 0 ;% 0.5 
rho_g = 0.8;

% sigma_u = (0.154/10000)*aa
 sigma_u = (0.154/100)*1
% sigma_u = (1.524/100)
sigma_g = (1.524/100)*2.0*0.8  %1.524

% Taylor Rule
a = 2.2;   % response to interest rate  
b = 0.5;   % response to GDP

A = [1 kappa; a*sigma b*sigma+1 ];
B = [beta 0; sigma 1];
% 
% eig(inv(B)*A)


A0 = [1 kappa; 0 1 ]; %É[Éçã‡óòÇ≈ÇÃçsóÒ


%% Step 1

% setting # of grids
N = 21; % # of grids of shocks
M = 3;  % # of grids of Gaussian Hermite



% N= max(N_u, N_g) 

u_range = mlow_u:(mup_u-mlow_u)/(N-1):mup_u;
g_range = mlow_g:(mup_g-mlow_g)/(N-1):mup_g;

[u,g] = meshgrid(u_range, g_range);

% setting policy functions
PF_pi = zeros(N,N);   %  inflation
PF_y = zeros(N,N);    % OUTPUT GAP
PF_i = zeros(N,N);    % INTEREST RATE 

%% Step 2

itr_num = 1;
  d = 0;  Tolerance = 0.0025; %1e-1; %

%   while d<1 
  while itr_num <= 40
 
     % save policy functions 
        old_PF_pi = PF_pi;
        old_PF_y = PF_y;
        old_PF_i = PF_i;
        
     % Derivation of Steady States
          pi_ss = interp2(u,g,old_PF_pi,0,0, 'cubic');
          y_ss  = interp2(u,g,old_PF_y, 0,0, 'cubic');
          i_ss = interp2(u,g,old_PF_i,0,0, 'cubic');
           
    
for j = 1:N       % # of grid g
    for k = 1:N  % # of grid u
                   
       %%  policy function of pi and y at t+1 period
            u_f = rho_u * u(j,k);   % j-> grid g, k->grid u
            g_f = rho_g * g(j,k);
            
            if forward_looking == 1  % perfect foresight 
              % linear interpolation in period t+1 
%                 pi_f = interp2(u,g,old_PF_pi,u_f,g_f, 'linear');
%                 y_f = interp2(u,g,old_PF_y,u_f,g_f, 'linear');
                
                pi_f = interp2(u,g,old_PF_pi,u_f,g_f, 'cubic')+pi_ss;
                y_f = interp2(u,g,old_PF_y,u_f,g_f, 'cubic')+y_ss;
                
            elseif forward_looking == 2    
               % stochastic 
                stochastic_forecast;  
            end
       
      %% calculation of endogenous variables in period t 
       %%  case of i >= -r_star
       
       
          i_temp = interp2(u,g,old_PF_i,u(j,k),g(j,k), 'linear') + i_ss;
           y = y_f -sigma*(i_temp-pi_f)+g(j,k) + y_ss;    % eq.3
           pi = beta*pi_f + kappa*y + u(j,k) + pi_ss;      % eq.2
           i =  a*pi + b*y +i_ss;               % eq.4 --> Taylor Rule
  
          
        
       %%  case of i < -r_star
        if (i < - r_star) & (const_ZLB==1)
             i = - r_star;                         % eq.4
            y = y_f -sigma*(i-pi_f)+g(j,k) + y_ss;         % eq.3
            pi = beta*pi_f + kappa*y + u(j,k) + pi_ss;      % eq.2
            
%              x = inv(A0)*B*[pi_f ; y_f] + inv(A0)*[0; -sigma]*r_star + inv(A0)*[u(j,k); g(j,k)];
%              pi = x(1); y=x(2);            
        end     
        
       %% update policy function
        ll = 0.1; ll_pi=0.1; ll_i = 0.2;
        PF_pi(j,k)= old_PF_pi(j,k)+ll*(pi-old_PF_pi(j,k));
        PF_y(j,k)= old_PF_y(j,k)+ll_pi*(y-old_PF_y(j,k));
%         PF_i(j,k)= i;
        PF_i(j,k)= old_PF_i(j,k)+ll_i*(i-old_PF_i(j,k)); 
       
      
         
    end
end   
      itr_num0 = num2str(itr_num);
      display([ itr_num0 'th-iteration']);
      itr_num = itr_num + 1;
   
   
%% Step 3  judgement of convergence

   Err_pi = old_PF_pi - PF_pi;
   Err_y = old_PF_y - PF_y;
   Err_i = old_PF_i - PF_i;
   
    err_pi = num2str(max(max(abs(Err_pi))));
    err_i = num2str(max(max(abs(Err_i))));
    err_y = num2str(max(max(abs(Err_y))));
    
   display( [ 'err_pi=' err_pi ',  err_i=' err_i ',  err_y='  err_y ] );
   
%    if (max(max(abs(Err_pi)))<Tolerance)&(max(max(abs(Err_y)))< Tolerance)&(max(max(abs(Err_i)))<Tolerance)
      if (max(max(abs(Err_pi)))< Tolerance)&(max(max(abs(Err_i)))<Tolerance)
         break;
%          d = 1;    
%    else
%        d = 0;
   end
   
end



%% Step 4Å@Making Graphs
%

figure(1)
surf(u,g,PF_pi)
% shading interp;
title('Policy Function of Inflation', 'Fontsize', 14);
   ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);       
      zlabel('Inflation', 'Fontsize', 14);
%       axis([-2 2 -2 2 -5 20]);
 
figure(2)
surf(u,g,PF_y)
% shading interp;
title('Policy Function of Output Gap', 'Fontsize', 14);
    ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);   
      zlabel('Output Gap', 'Fontsize', 14);
      
figure(3)
surf(u,g,PF_i)
% shading interp;
title('Policy Function of Interest Rate', 'Fontsize', 14);
    ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);    
    zlabel('Interest Rate', 'Fontsize', 14);
         
figure(4)
subplot(2,2,1)
surf(u,g,Err_pi)
% shading interp;
title('Convergence of Policy Function of Inflation', 'Fontsize', 12);
   ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);   
      zlabel('Err', 'Fontsize', 14);
%       axis([-0.5 0.5 -0.5 0.5 -0.1 0.1]);

subplot(2,2,2)
surf(u,g,Err_y)
% shading interp;
title('Convergence of Policy Function of Output Gap', 'Fontsize', 12);
  ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);    
      zlabel('Err', 'Fontsize', 14);

subplot(2,2,3)      
surf(u,g,Err_i)
% shading interp;
title('Convergence of Policy Function of Interest Rate', 'Fontsize', 12);
     ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);   
      zlabel('Err', 'Fontsize', 14);   


%% Impulse Response Functions
nsim = 100;    % # of Monte Calro
hh = 20;       % horizons of response  

%% positive markup shock
u_sim1 = 0.05; % markup shock
g_sim1 = 0;    % real rate shock
pi_sim1 = 0; 
y_sim1 = 0; %  
i_sim1 = 0; %

pi0 = interp2(u, g, PF_pi, 0, 0, 'linear');
y0  = interp2(u, g, PF_y, 0, 0, 'linear');
i0  = interp2(u, g, PF_i, 0, 0, 'linear');

for j = 1:nsim
  %   
  pi_sim1(j) = interp2(u, g, PF_pi, u_sim1(j), g_sim1(j), 'linear');
  y_sim1(j) = interp2(u, g, PF_y, u_sim1(j), g_sim1(j), 'linear');
  i_sim1(j) = interp2(u, g, PF_i, u_sim1(j), g_sim1(j), 'linear');
  
  if j < nsim
    u_sim1(j+1) = rho_u * u_sim1(j);
    g_sim1(j+1) = rho_g * g_sim1(j);
  end
end

figure(5)
subplot(3,1,3)
plot(0:hh, [y0,y_sim1(1:hh)]); title('y')
subplot(3,1,2)
plot(0:hh, [pi0,pi_sim1(1:hh)]); title('pi')
subplot(3,1,1)
plot(0:hh, [i0,i_sim1(1:hh)]); title('i')
%subplot(2,2,4)
%plot(0:50, [0,u_sim1(1:50)]); title('u')

%% negative markup shock

u_sim2 = -0.05;  % markup shock
g_sim2 = 0;      % real rate shock
pi_sim2 = 0; 
y_sim2 = 0; %  
i_sim2 = 0; %

pi0 = interp2(u, g, PF_pi, 0, 0, 'linear');
y0  = interp2(u, g, PF_y, 0, 0, 'linear');
i0  = interp2(u, g, PF_i, 0, 0, 'linear');

for j = 1:nsim
  pi_sim2(j) = interp2(u, g, PF_pi, u_sim2(j), g_sim2(j), 'linear');
  y_sim2(j) = interp2(u, g, PF_y, u_sim2(j), g_sim2(j), 'linear');
  i_sim2(j) = interp2(u, g, PF_i, u_sim2(j), g_sim2(j), 'linear');
%   if j > 1
%     i_sim2(j-1) = max(sigma^(-1)*(y_sim2(j) - y_sim2(j-1) + sigma*pi_sim2(j)+g_sim2(j-1)), -r_star);
%   end
  if j < nsim
    u_sim2(j+1) = rho_u * u_sim2(j);
    g_sim2(j+1) = rho_g * g_sim2(j);
  end
end

figure(6)
   y_sim = [y0 y_sim1(1:hh);    y0 y_sim2(1:hh) ];
   pi_sim = [pi0 pi_sim1(1:hh); pi0 pi_sim2(1:hh) ];
   i_sim = [i0 i_sim1(1:hh);    i0 i_sim2(1:hh) ];
subplot(3,1,3)   
plot(0:hh, y_sim,'LineWidth',2 ); title('y'); 
legend('positive markup shock', 'negative markup shock'); 
subplot(3,1,2)
plot(0:hh, pi_sim,'LineWidth',2); title('pi')
legend('positive markup shock', 'negative markup shock'); 
subplot(3,1,1)
plot(0:hh, i_sim,'LineWidth', 2  ); title('i')
legend('positive markup shock', 'negative markup shock'); 


%%
if (flag_save_policy_func == 1)&(forward_looking == 1) &(const_ZLB==1)
     PF_y_Taylor_pf = PF_y;
     PF_i_Taylor_pf = PF_i;
     PF_pi_Taylor_pf = PF_pi;
     Taylor_pf_u = u;
     Taylor_pf_g = g;
     file ='Taylor_pf_ZLB_no_switch.mat';   
     save('./output/Taylor_pf_ZLB_no_switch.mat', 'PF_y_Taylor_pf','PF_i_Taylor_pf','PF_pi_Taylor_pf','Taylor_pf_u','Taylor_pf_g');

elseif (flag_save_policy_func == 1)&(forward_looking == 2) &(const_ZLB==1)
     PF_y_Taylor_stc = PF_y;
     PF_i_Taylor_stc = PF_i;
     PF_pi_Taylor_stc = PF_pi;
     Taylor_stc_u = u;
     Taylor_stc_g = g;
     file ='Taylor_stc_ZLB_no_switch.mat';   
     save('./output/Taylor_stc_ZLB_no_switch.mat', 'PF_y_Taylor_stc','PF_i_Taylor_stc','PF_pi_Taylor_stc','Taylor_stc_u','Taylor_stc_g');
elseif    (flag_save_policy_func == 1)&(forward_looking == 1) 
     PF_y_Taylor_pf = PF_y;
     PF_i_Taylor_pf = PF_i;
     PF_pi_Taylor_pf = PF_pi;
     Taylor_pf_u = u;
     Taylor_pf_g = g;
     file ='Taylor_pf_no_switch.mat';   
     save('./output/Taylor_pf_no_switch.mat', 'PF_y_Taylor_pf','PF_i_Taylor_pf','PF_pi_Taylor_pf','Taylor_pf_u','Taylor_pf_g');
else
     PF_y_Taylor_stc = PF_y;
     PF_i_Taylor_stc = PF_i;
     PF_pi_Taylor_stc = PF_pi;
     Taylor_stc_u = u;
     Taylor_stc_g = g;
     file ='Taylor_stc_no_switch.mat';   
     save('./output/Taylor_stc_no_switch.mat', 'PF_y_Taylor_stc','PF_i_Taylor_stc','PF_pi_Taylor_stc','Taylor_stc_u','Taylor_stc_g');
end     
    
     