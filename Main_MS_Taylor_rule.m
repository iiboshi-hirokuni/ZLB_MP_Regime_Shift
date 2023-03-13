%%
%%    By Hirokuni Iiboshi (2016)
%% @"Monetary Policy Regime Shifts under the Zero Lower Bound:
%%    An application of a stochastic rational expectations equilibrium to a Markov switching DSGE model" 
%%    Economic Modelling, Vol. 52, p186-205  
%%
%%

addpath('functions/');

%%  save files of policy functions
flag_save_policy_func = 1; % 1-> on,  0-> off

%%  Option of ZLB
const_ZLB = 1;   %  1-> ZLB,  0-> non ZLB

forward_looking = 2 % 1: perfect foresight, 2: stochastic  


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
rho_u = 0.0 ; 
rho_g = 0.8;

sigma_u = (0.1524/100)*1
sigma_g = (1.524/100)*2.0  %1.524

% Taylor Rule
% regime 0 (Aggressive)
a0 = 2.2;   % response to inflation  
b0 = 0.5;   % response to GDP

% regiem 1 (Passive)
a1 = 0.8;   % response to inflation
b1 = 0.15;   % response to GDP

% Transtion Prob
p00 = 0.75;    % Prob( St = 0 | St-1 = 0) (Aggressive)
p11 = 0.75;    % Prob( St = 1 | St-1 = 1) (Passive)

P = [p00 1-p11; ...
     1-p00 p11];


%% Step 1

% setting # of grids
N = 21; % # of grids of shocks
M = 3;  % # of grids of Gaussian Hermite


% N= max(N_u, N_g) 

u_range = mlow_u:(mup_u-mlow_u)/(N-1):mup_u;
g_range = mlow_g:(mup_g-mlow_g)/(N-1):mup_g;

[u,g] = meshgrid(u_range, g_range);

% setting policy functions
% Regiem 0
 PF0_pi = zeros(N,N);   % inflation
 PF0_y = zeros(N,N);    % OUTPUT GAP
 PF0_i = zeros(N,N);    % INTEREST RATE
% Regime 1
 PF1_pi = zeros(N,N);   % inflation
 PF1_y = zeros(N,N);    % OUTPUT GAP
 PF1_i = zeros(N,N);    % INTEREST RATE

 
%% Step 2

itr_num = 1;
  d = 0;  Tolerance = 0.0025; %1e-1; %

%   while d<1 
  while itr_num <= 40
 
     % save policy functions 
     % regime 0
        old_PF0_pi = PF0_pi;
        old_PF0_y = PF0_y;
        old_PF0_i = PF0_i;
      % regime 1
        old_PF1_pi = PF1_pi;
        old_PF1_y = PF1_y;
        old_PF1_i = PF1_i;
        
      % Derivation of Steady States
          pi_ss0 = interp2(u,g,old_PF0_pi,0,0, 'cubic');
          y_ss0  = interp2(u,g,old_PF0_y, 0,0, 'cubic');
          i_ss0 = interp2(u,g,old_PF0_i,0,0, 'cubic');
          pi_ss1 = interp2(u,g,old_PF1_pi,0,0, 'cubic');
          y_ss1  = interp2(u,g,old_PF1_y, 0,0, 'cubic');
          i_ss1 = interp2(u,g,old_PF1_i,0,0, 'cubic');
        
for j = 1:N      % # of grid g
    for k = 1:N  % # of grid u
        
       %%  policy function of pi and y at t+1 period
            u_f = rho_u * u(j,k);   % j-> grid g, k->grid u
            g_f = rho_g * g(j,k);
            
            if forward_looking == 1  % perfect foresight 
              % linear interpolation in period t+1            
               pi_f0 = interp2(u,g,old_PF0_pi,u_f,g_f, 'cubic') + pi_ss0;
               pi_f1 = interp2(u,g,old_PF1_pi,u_f,g_f, 'cubic') + pi_ss1 ;                       
               
               y_f0 = interp2(u,g,old_PF0_y,u_f,g_f, 'cubic') + y_ss0;
               y_f1 = interp2(u,g,old_PF1_y,u_f,g_f, 'cubic') + y_ss1;               
                 
            % regime 0
                pi0_f = P(1,1)*pi_f0 +  P(2,1)*pi_f1;   
                y0_f = P(1,1)*y_f0 +  P(2,1)*y_f1; 
            % regime 1
                pi1_f = P(1,2)*pi_f0 +  P(2,2)*pi_f1;   
                y1_f = P(1,2)*y_f0 +  P(2,2)*y_f1; 
                
               
            elseif forward_looking == 2    
               % stochastic 
                 stochastic_forecast_MS;  
            end
       
       %% calculation of endogenous variables in period t 
       %%  case of i >= -r_star      
       
       %% regime 0           
           i0_temp = interp2(u,g,old_PF0_i,u(j,k),g(j,k), 'linear')+ i_ss0;
           y0 = y0_f -sigma*(i0_temp-pi0_f)+g(j,k) + y_ss0 ;    % eq.3
           pi0 = beta*pi0_f + kappa*y0 + u(j,k)+ pi_ss0  ;      % eq.2
           i0 =   a0*pi0 + b0*y0 + i_ss0 ;               % eq.4 --> Taylor Rule
           
           %%  case of i < -r_star
        if (i0 < - r_star) & (const_ZLB==1)
            i0 = - r_star;                         % eq.4
            y0 = y0_f -sigma*(i0-pi0_f)+g(j,k) + y_ss0 ;          % eq.3
            pi0 = beta*pi0_f + kappa*y0 + u(j,k)+ pi_ss0 ;       % eq.2       
        end     
        
  
        %% regime 1           
           i1_temp = interp2(u,g,old_PF1_i,u(j,k),g(j,k), 'linear')+ i_ss1;
           y1 = y1_f -sigma*(i1_temp-pi1_f)+g(j,k)+ y_ss1;    % eq.3
           pi1 = beta*pi1_f + kappa*y1 + u(j,k)+ pi_ss1;      % eq.2
           i1 =  a1*pi1 + b1*y1 + i_ss1;              % eq.4 --> Taylor Rule
        
        %%  case of i < -r_star
        if (i1 < - r_star) & (const_ZLB==1)
            i1 = - r_star;                         % eq.4
            y1 = y1_f -sigma*(i1-pi1_f)+g(j,k)+ y_ss1;         % eq.3
            pi1 = beta*pi1_f + kappa*y1 + u(j,k)+ pi_ss1;      % eq.2       
        end     
        
        
       %% update policy function
        ll = 0.1; ll_pi=0.1; ll_i = 0.2;
        %% regime 0   
        PF0_pi(j,k)= old_PF0_pi(j,k)+ll*(pi0-old_PF0_pi(j,k));
        PF0_y(j,k)= old_PF0_y(j,k)+ll_pi*(y0-old_PF0_y(j,k));
        PF0_i(j,k)= old_PF0_i(j,k)+ll_i*(i0-old_PF0_i(j,k)); 
        %% regime 1   
        PF1_pi(j,k)= old_PF1_pi(j,k)+ll*(pi1-old_PF1_pi(j,k));
        PF1_y(j,k)= old_PF1_y(j,k)+ll_pi*(y1-old_PF1_y(j,k));
        PF1_i(j,k)= old_PF1_i(j,k)+ll_i*(i1-old_PF1_i(j,k)); 
         
    end
end   
      itr_num0 = num2str(itr_num);
      display([ itr_num0 'th-iteration']);
      itr_num = itr_num + 1;
   
   
%% Step 3  judgement of convergence

   Err0_pi = old_PF0_pi - PF0_pi;    
   Err0_y = old_PF0_y - PF0_y;     
   Err0_i = old_PF0_i - PF0_i;
   
   Err1_pi = old_PF1_pi - PF1_pi;
   Err1_y = old_PF1_y - PF1_y;
   Err1_i = old_PF1_i - PF1_i;
   
    err0_pi = num2str(max(max(abs(Err0_pi))));
    err0_i = num2str(max(max(abs(Err0_i))));
    err0_y = num2str(max(max(abs(Err0_y))));
    
    err1_pi = num2str(max(max(abs(Err1_pi))));
    err1_i = num2str(max(max(abs(Err1_i))));
    err1_y = num2str(max(max(abs(Err1_y))));
    
   display( [ 'err0_pi=' err0_pi ',  err0_i=' err0_i ',  err0_y='  err0_y ] );
   display( [ 'err1_pi=' err1_pi ',  err1_i=' err1_i ',  err1_y='  err1_y ] );
   
%    if (max(max(abs(Err_pi)))<Tolerance)&(max(max(abs(Err_y)))< Tolerance)&(max(max(abs(Err_i)))<Tolerance)
      if (max(max(abs(Err0_pi)))< Tolerance)&(max(max(abs(Err0_i)))<Tolerance)
         break;
%          d = 1;    
%    else
%        d = 0;
   end
   
end



%% Step 4.  Making Graphs
%

figure(1)
subplot(1,2,1)
surf(u,g,PF0_pi)
% shading interp;
title('Policy Function of Inflation: Regime 0', 'Fontsize', 14);
      ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);      
      zlabel('Inflation', 'Fontsize', 14);
subplot(1,2,2)
surf(u,g,PF1_pi)
% shading interp;
title('Policy Function of Inflation: Regime 1', 'Fontsize', 14);
      ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);        
      zlabel('Inflation', 'Fontsize', 14);
 
figure(2)
subplot(1,2,1)
surf(u,g,PF0_y)
% shading interp;
title('Policy Function of Output Gap: Regime 0', 'Fontsize', 14);
     ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);   
      zlabel('Output Gap', 'Fontsize', 14);
subplot(1,2,2)
surf(u,g,PF1_y)
% shading interp;
title('Policy Function of Output Gap: Regime 1', 'Fontsize', 14);
     ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);        
     zlabel('Output Gap', 'Fontsize', 14); 
      
figure(3)
subplot(1,2,1)
surf(u,g,PF0_i)
% shading interp;
title('Policy Function of Interest Rate: Regime 0', 'Fontsize', 14);
    ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);       
    zlabel('Interest Rate', 'Fontsize', 14);
subplot(1,2,2)
surf(u,g,PF1_i)
% shading interp;
title('Policy Function of Interest Rate: Regime 1', 'Fontsize', 14);
    ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);         
      zlabel('Inflation', 'Fontsize', 14);   
      
figure(4)
subplot(2,2,1)
surf(u,g,Err0_pi)
% shading interp;
title('Convergence of Policy Function of Inflation', 'Fontsize', 12);
    ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);   
      zlabel('Err', 'Fontsize', 14);
%       axis([-0.5 0.5 -0.5 0.5 -0.1 0.1]);

subplot(2,2,2)
surf(u,g,Err0_y)
% shading interp;
title('Convergence of Policy Function of Output Gap', 'Fontsize', 12);
   ylabel('g_t(real rate shock)', 'Fontsize', 12);
      xlabel('u_t(markup shock)', 'Fontsize', 12);   
      zlabel('Err', 'Fontsize', 14);

subplot(2,2,3)      
surf(u,g,Err0_i)
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

 pi_sim0 = interp2(u, g, PF0_pi, 0, 0, 'linear');
  y_sim0 = interp2(u, g, PF0_y, 0, 0, 'linear');
  i_sim0 = interp2(u, g, PF0_i, 0, 0, 'linear');

for j = 1:nsim
  pi_sim1(j) = interp2(u, g, PF0_pi, u_sim1(j), g_sim1(j), 'linear');
  y_sim1(j) = interp2(u, g, PF0_y, u_sim1(j), g_sim1(j), 'linear');
  i_sim1(j) = interp2(u, g, PF0_i, u_sim1(j), g_sim1(j), 'linear');
  
  if j < nsim
    u_sim1(j+1) = rho_u * u_sim1(j);
    g_sim1(j+1) = rho_g * g_sim1(j);
  end
end


%% negative markup shock

u_sim2 = -0.05;  % markup shock
g_sim2 = 0;      % real rate shock
pi_sim2 = 0; 
y_sim2 = 0; %  
i_sim2 = 0; %

for j = 1:nsim
  pi_sim2(j) = interp2(u, g, PF0_pi, u_sim2(j), g_sim2(j), 'linear');
  y_sim2(j) = interp2(u, g, PF0_y, u_sim2(j), g_sim2(j), 'linear');
  i_sim2(j) = interp2(u, g, PF0_i, u_sim2(j), g_sim2(j), 'linear');
%   if j > 1
%     i_sim2(j-1) = max(sigma^(-1)*(y_sim2(j) - y_sim2(j-1) + sigma*pi_sim2(j)+g_sim2(j-1)), -r_star);
%   end
  if j < nsim
    u_sim2(j+1) = rho_u * u_sim2(j);
    g_sim2(j+1) = rho_g * g_sim2(j);
  end
end

figure(5)
   y_sim = [y_sim0 y_sim1(1:hh);    y_sim0 y_sim2(1:hh) ];
   pi_sim = [pi_sim0 pi_sim1(1:hh); pi_sim0 pi_sim2(1:hh) ];
   i_sim = [i_sim0 i_sim1(1:hh);    i_sim0 i_sim2(1:hh) ];
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
     PF0_y_Taylor_pf = PF0_y;
     PF0_i_Taylor_pf = PF0_i;
     PF0_pi_Taylor_pf = PF0_pi;
     PF1_y_Taylor_pf = PF1_y;
     PF1_i_Taylor_pf = PF1_i;
     PF1_pi_Taylor_pf = PF1_pi;     
     Taylor_pf_u = u;
     Taylor_pf_g = g;
     file ='Taylor_pf_ZLB.mat';   
     save('./output/Taylor_pf_ZLB.mat', 'PF0_y_Taylor_pf','PF0_i_Taylor_pf', 'PF0_pi_Taylor_pf', ...
     'PF1_y_Taylor_pf','PF1_i_Taylor_pf', 'PF1_pi_Taylor_pf', 'Taylor_pf_u','Taylor_pf_g');
elseif (flag_save_policy_func == 1)&(forward_looking == 2) &(const_ZLB==1)
     PF0_y_Taylor_stc = PF0_y;
     PF0_i_Taylor_stc = PF0_i;
     PF0_pi_Taylor_stc = PF0_pi;
     PF1_y_Taylor_stc = PF1_y;
     PF1_i_Taylor_stc = PF1_i;
     PF1_pi_Taylor_stc = PF1_pi;     
     Taylor_stc_u = u;
     Taylor_stc_g = g;
     file ='Taylor_stc_ZLB.mat';   
     save('./output/Taylor_stc_ZLB.mat', 'PF0_y_Taylor_stc','PF0_i_Taylor_stc','PF0_pi_Taylor_stc',...
         'PF1_y_Taylor_stc','PF1_i_Taylor_stc','PF1_pi_Taylor_stc', 'Taylor_stc_u','Taylor_stc_g');
elseif    (flag_save_policy_func == 1)&(forward_looking == 1) 
     PF0_y_Taylor_pf = PF0_y;
     PF0_i_Taylor_pf = PF0_i;
     PF0_pi_Taylor_pf = PF0_pi;
     PF1_y_Taylor_pf = PF1_y;
     PF1_i_Taylor_pf = PF1_i;
     PF1_pi_Taylor_pf = PF1_pi;
     Taylor_pf_u = u;
     Taylor_pf_g = g;
     file ='Taylor_pf.mat';   
     save('./output/Taylor_pf.mat', 'PF0_y_Taylor_pf','PF0_i_Taylor_pf','PF0_pi_Taylor_pf',...
         'PF1_y_Taylor_pf','PF1_i_Taylor_pf','PF1_pi_Taylor_pf','Taylor_pf_u','Taylor_pf_g');
else
     PF0_y_Taylor_stc = PF0_y;
     PF0_i_Taylor_stc = PF0_i;
     PF0_pi_Taylor_stc = PF0_pi;
     PF1_y_Taylor_stc = PF1_y;
     PF1_i_Taylor_stc = PF1_i;
     PF1_pi_Taylor_stc = PF1_pi;     
     Taylor_stc_u = u;
     Taylor_stc_g = g;
     file ='Taylor_stc.mat';   
     save('./output/Taylor_stc.mat', 'PF0_y_Taylor_stc','PF0_i_Taylor_stc','PF0_pi_Taylor_stc',...
           'PF1_y_Taylor_stc','PF1_i_Taylor_stc','PF1_pi_Taylor_stc', 'Taylor_stc_u','Taylor_stc_g');
end     
    
     