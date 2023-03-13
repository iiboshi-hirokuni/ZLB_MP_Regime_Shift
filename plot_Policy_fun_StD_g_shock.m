%%
%%    By Hirokuni Iiboshi (2016)
%% @"Monetary Policy Regime Shifts under the Zero Lower Bound:
%%    An application of a stochastic rational expectations equilibrium to a Markov switching DSGE model" 
%%    Economic Modelling, Vol. 52, p186-205  
%%
%%

% setting # of grids
N = 21; % # of grids of shocks

% setting ranges of structural shocks
mlow_u = -0.1; mup_u = 0.1;  % markup shock
N_u = round((mup_u - mlow_u)/0.005);

mlow_g = -0.4; mup_g = 0.4; % real rate shock
N_g = round((mup_g - mlow_g)/0.002);
 
 low_g =-0.1;
 up_g  = 0.1;
% low_g =mlow_g;
% up_g  =mup_g;


 g_sim = mlow_g:(mup_g-mlow_g)/(N-1):mup_g;
 u_sim = -0.00;
 
%======================================  
%
%=====================================
regime = 0;  % Aggressive->0 passive->1

  
  model = 'Taylor_ZLB'
%  model = 'Taylor'

switch model
 
 case 'Taylor_ZLB' 
      load('./output_2/Taylor_stc_ZLB.mat', 'PF0_y_Taylor_stc');
     load('./output_2/Taylor_stc_ZLB.mat', 'PF0_i_Taylor_stc');
     load('./output_2/Taylor_stc_ZLB.mat', 'PF0_pi_Taylor_stc');
     load('./output_2/Taylor_stc_ZLB.mat', 'PF1_y_Taylor_stc');
     load('./output_2/Taylor_stc_ZLB.mat', 'PF1_i_Taylor_stc');
     load('./output_2/Taylor_stc_ZLB.mat', 'PF1_pi_Taylor_stc');     
     load('./output_2/Taylor_stc_ZLB.mat', 'Taylor_stc_u');
     load('./output_2/Taylor_stc_ZLB.mat', 'Taylor_stc_g');
     
%      load('./output/Taylor_pf_ZLB.mat', 'Taylor_pf_u');
%      load('./output/Taylor_pf_ZLB.mat', 'Taylor_pf_g');

      if regime ==1
      PF_y_pf = PF1_y_Taylor_stc;
      PF_i_pf = PF1_i_Taylor_stc;
      PF_pi_pf = PF1_pi_Taylor_stc;
      else
     PF_y_pf = PF0_y_Taylor_stc;
     PF_i_pf = PF0_i_Taylor_stc;
     PF_pi_pf = PF0_pi_Taylor_stc;
      end
      
     u_pf = Taylor_stc_u;
     g_pf = Taylor_stc_g;
     
load('./output/Taylor_stc_ZLB.mat', 'PF0_y_Taylor_stc');
load('./output/Taylor_stc_ZLB.mat', 'PF0_i_Taylor_stc');
load('./output/Taylor_stc_ZLB.mat', 'PF0_pi_Taylor_stc');
load('./output/Taylor_stc_ZLB.mat', 'PF1_y_Taylor_stc');
load('./output/Taylor_stc_ZLB.mat', 'PF1_i_Taylor_stc');
load('./output/Taylor_stc_ZLB.mat', 'PF1_pi_Taylor_stc');

load('./output/Taylor_stc_ZLB.mat', 'Taylor_stc_u');
load('./output/Taylor_stc_ZLB.mat', 'Taylor_stc_g');

  if regime ==1 
      PF_y_stc = PF1_y_Taylor_stc;
      PF_i_stc = PF1_i_Taylor_stc;
      PF_pi_stc = PF1_pi_Taylor_stc;
  else    
     PF_y_stc = PF0_y_Taylor_stc;
     PF_i_stc = PF0_i_Taylor_stc;
     PF_pi_stc = PF0_pi_Taylor_stc;
  end    
     u_stc = Taylor_stc_u;
     g_stc = Taylor_stc_g;
        
 

case 'Taylor' 
     load('./output/Taylor_pf.mat', 'PF_y_Taylor_pf');
     load('./output/Taylor_pf.mat', 'PF_i_Taylor_pf');
     load('./output/Taylor_pf.mat', 'PF_pi_Taylor_pf');
     load('./output/Taylor_pf.mat', 'Taylor_pf_u');
     load('./output/Taylor_pf.mat', 'Taylor_pf_g');

     PF_y_pf = PF_y_Taylor_pf;
     PF_i_pf = PF_i_Taylor_pf;
     PF_pi_pf = PF_pi_Taylor_pf;
     u_pf = Taylor_pf_u;
     g_pf = Taylor_pf_g;
     
load('./output/Taylor_stc.mat', 'PF_y_Taylor_stc');
load('./output/Taylor_stc.mat', 'PF_i_Taylor_stc');
load('./output/Taylor_stc.mat', 'PF_pi_Taylor_stc');
load('./output/Taylor_stc.mat', 'Taylor_stc_u');
load('./output/Taylor_stc.mat', 'Taylor_stc_g');

     PF_y_stc = PF_y_Taylor_stc;
     PF_i_stc = PF_i_Taylor_stc;
     PF_pi_stc = PF_pi_Taylor_stc;
     u_stc = Taylor_stc_u;
     g_stc = Taylor_stc_g;
  
     
end

 
%   save2_PF_pi = save_PF_pi;
%   save2_PF_y = save_PF_y;
%   save2_PF_i = save_PF_i;
%   save2_u = save_u; save2_g = save_g;

%  save_PF_pi = PF_pi;
%  save_PF_y = PF_y;
%  save_PF_i = PF_i;
%  save_u = u; save_g = g;
 
 y=zeros(N,1); i=zeros(N,1); pi=zeros(N,1);
 y2=zeros(N,1); i2=zeros(N,1); pi2=zeros(N,1);
 
%   u_sim = -0.002;
%   u_sim = -0.00;

 pi_0 = interp2(u_stc, g_stc, PF_pi_stc, 0, 0, 'cubic');
 i_0 = interp2(u_stc, g_stc, PF_i_stc, 0, 0, 'cubic');
 y_0 = interp2(u_stc, g_stc, PF_y_stc, 0, 0, 'cubic');
 pi_ss = interp2(u_pf, g_pf, PF_pi_pf, 0, 0, 'cubic');
 i_ss = interp2(u_pf, g_pf, PF_i_pf, 0, 0, 'cubic');
 y_ss = interp2(u_pf, g_pf, PF_y_pf, 0, 0, 'cubic');
  
  for j = 1:N
     
     pi2(j) = interp2(u_pf, g_pf, PF_pi_pf, u_sim, g_sim(j), 'cubic')-pi_ss;
     y2(j)  = interp2(u_pf, g_pf, PF_y_pf, u_sim, g_sim(j), 'cubic')-y_ss;
     i2(j)  = interp2(u_pf, g_pf, PF_i_pf, u_sim, g_sim(j), 'cubic');
     
       pi(j) = interp2(u_stc, g_stc, PF_pi_stc, u_sim, g_sim(j), 'cubic')-pi_0;
       y(j)  = interp2(u_stc, g_stc, PF_y_stc, u_sim, g_sim(j), 'cubic')-y_0;
       i(j)  = interp2(u_stc, g_stc, PF_i_stc, u_sim, g_sim(j), 'cubic');
     
end

  

figure(10)
  
  subplot(3,1,1)
     hold on
       plot(100*g_sim,100*[y'],'b','LineWidth',2)
        plot(100*g_sim,100*[y2'],'r--','LineWidth',2)
%          ylim([-0.1, 0.02]);
       xlim(100*[low_g,up_g])
       title('Output','Fontsize', 14); 
      xlabel('g_t (Real Interest Rate Shock)', 'Fontsize', 12);

      legend('20 % reduction of StD', 'Original')
      
%        legend('stochastic*0.5', 'parfect forsight')
%         legend( 'stochastic of NK', 'parfect forsight of NK' )
     hold off  
  subplot(3,1,2)
       hold on
       plot(100*g_sim,100*[pi'],'b','LineWidth',2)
        plot(100*g_sim,100*[pi2'],'r--','LineWidth',2)
%        plot(g_sim,[pi'; pi2'],'LineWidth',2)
         %ylim([-0.01, 0.005]);
        xlim(100*[low_g,up_g])
%        xlim([-0.08,0.02])
       title('Inflation', 'Fontsize', 14); 
      xlabel('g_t (Real Interest Rate Shock)', 'Fontsize', 12);
       legend('20 % reduction of StD', 'Original')
       hold off;
  subplot(3,1,3)
      hold on
       plot(100*g_sim,100*[i'],'b','LineWidth',2)
        plot(100*g_sim,100*[i2'],'r--','LineWidth',2)
%        plot(g_sim,[i'; i2'],'LineWidth',2)
%         ylim([-0.01, 0.02]);
         xlim(100*[low_g,up_g])
%          xlim([-0.08,0.02])
       title('Nominal Interest Rate', 'Fontsize', 14); 
       xlabel('g_t (Real Interest Rate Shock)', 'Fontsize', 12);
       legend('20 % reduction of StD', 'Original')
       hold off
