%%
%%    By Hirokuni Iiboshi (2016)
%% Å@"Monetary Policy Regime Shifts under the Zero Lower Bound:
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
%  
 
 set =1;
if( set==1)
     low_g = -0.1; %-0.06;
     up_g  = 0.1;  % -0.0;
else 
     low_g =mlow_g;
     up_g  =mup_g;
end     


 g_sim = mlow_g:(mup_g-mlow_g)/(N-1):mup_g;
 
  %  model = 'Taylor_pf'
    model = 'Taylor_stc'

switch model
 
 case 'Taylor_pf' 
     load('./output/Taylor_pf_ZLB.mat', 'PF0_y_Taylor_pf');
     load('./output/Taylor_pf_ZLB.mat', 'PF0_i_Taylor_pf');
     load('./output/Taylor_pf_ZLB.mat', 'PF0_pi_Taylor_pf');
    load('./output/Taylor_pf_ZLB.mat', 'PF1_y_Taylor_pf');
     load('./output/Taylor_pf_ZLB.mat', 'PF1_i_Taylor_pf');
     load('./output/Taylor_pf_ZLB.mat', 'PF1_pi_Taylor_pf');
     
     load('./output/Taylor_pf_ZLB.mat', 'Taylor_pf_u');
     load('./output/Taylor_pf_ZLB.mat', 'Taylor_pf_g');

     PF_y_ZLB = PF0_y_Taylor_pf;
     PF_i_ZLB = PF0_i_Taylor_pf;
     PF_pi_ZLB = PF0_pi_Taylor_pf;
     u_ZLB = Taylor_pf_u;
     g_ZLB = Taylor_pf_g;
   

     PF_y_no_ZLB = PF1_y_Taylor_pf;
     PF_i_no_ZLB = PF1_i_Taylor_pf;
     PF_pi_no_ZLB = PF1_pi_Taylor_pf;
     u_no_ZLB = Taylor_pf_u;
     g_no_ZLB = Taylor_pf_g;
 
case 'Taylor_stc' 
     load('./output/Taylor_stc_ZLB.mat', 'PF0_y_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF0_i_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF0_pi_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF1_y_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF1_i_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF1_pi_Taylor_stc');     
     load('./output/Taylor_stc_ZLB.mat', 'Taylor_stc_u');
     load('./output/Taylor_stc_ZLB.mat', 'Taylor_stc_g');

     PF_y_ZLB = PF0_y_Taylor_stc;
     PF_i_ZLB = PF0_i_Taylor_stc;
     PF_pi_ZLB = PF0_pi_Taylor_stc;
     u_ZLB = Taylor_stc_u;
     g_ZLB = Taylor_stc_g;
  
     PF_y_no_ZLB = PF1_y_Taylor_stc;
     PF_i_no_ZLB = PF1_i_Taylor_stc;
     PF_pi_no_ZLB = PF1_pi_Taylor_stc;
      u_no_ZLB = Taylor_stc_u;
      g_no_ZLB = Taylor_stc_g;
        
     
     
end

% mlow = -0.1; mup = 0.0;

 
 y=zeros(N,1); i=zeros(N,1); pi=zeros(N,1);
 y2=zeros(N,1); i2=zeros(N,1); pi2=zeros(N,1);
 
%   u_sim = -0.002;
  u_sim = -0.00;
  
  for j = 1:N
     
     pi2(j) = interp2(u_no_ZLB, g_no_ZLB, PF_pi_no_ZLB, u_sim, g_sim(j), 'cubic');
     y2(j)  = interp2(u_no_ZLB, g_no_ZLB, PF_y_no_ZLB, u_sim, g_sim(j), 'cubic');
     i2(j)  = interp2(u_no_ZLB, g_no_ZLB, PF_i_no_ZLB, u_sim, g_sim(j), 'cubic');
     
       pi(j) = interp2(u_ZLB, g_ZLB, PF_pi_ZLB, u_sim, g_sim(j), 'cubic');
       y(j)  = interp2(u_ZLB, g_ZLB, PF_y_ZLB, u_sim, g_sim(j), 'cubic');
       i(j)  = interp2(u_ZLB, g_ZLB, PF_i_ZLB, u_sim, g_sim(j), 'cubic');
     
end

  

figure
  
  subplot(3,1,1)
     hold on
       plot(100*g_sim,100*[y'],'b','LineWidth',2)
        plot(100*g_sim,100*[y2'],'r--','LineWidth',2)
        ylim([-10.0, 10.0]);
       xlim(100*[low_g,up_g])
       title('Output', 'fontsize', 15); 
      %xlabel('g_t (Real Interest Rate Shock)', 'Fontsize', 12);
      xlabel('g_t (Aggregate Demand Shock)', 'Fontsize', 12);

      switch model
       case 'Taylor_pf' %|'Taylor_ZLB'
         legend('Active Regime', 'Passive Regime')  
       case 'Taylor_stc'
         legend('Active Regime', 'Passive Regime')           
       case 'Opt_pf' 
         legend('ZLB of Optimal in parfect forsight', 'no ZLB')
       case 'Opt_stc' 
         legend('ZLB of Optimal in stochastic', 'no ZLB')
      end  
%        legend('stochastic*0.5', 'parfect forsight')
%         legend( 'stochastic of NK', 'parfect forsight of NK' )
     hold off  
  subplot(3,1,2)
       hold on
       plot(100*g_sim,100*[pi'],'b','LineWidth',2)
        plot(100*g_sim,100*[pi2'],'r--','LineWidth',2)
%        plot(g_sim,[pi'; pi2'],'LineWidth',2)
%          ylim([-0.01, 0.005]);
        xlim(100*[low_g,up_g])
%        xlim([-0.08,0.02])
      title('Inflation', 'fontsize', 15); 
        %xlabel('g_t (Real Interest Rate Shock)', 'Fontsize', 12);
      xlabel('g_t (Aggregate Demand Shock)', 'Fontsize', 12);
      
        switch model
        case 'Taylor_pf' %|'Taylor_ZLB'
         legend('Active Regime', 'Passive Regime')  
       case 'Taylor_stc'
         legend('Active Regime', 'Passive Regime')         
       case 'Opt_pf' 
         legend('ZLB of Optimal in parfect forsight', 'no ZLB')
       case 'Opt_stc' 
         legend('ZLB of Optimal in stochastic', 'no ZLB')
      end  
       hold off;
  subplot(3,1,3)
      hold on
       plot(100*g_sim,100*[i'],'b','LineWidth',2)
        plot(100*g_sim,100*[i2'],'r--','LineWidth',2)
%        plot(g_sim,[i'; i2'],'LineWidth',2)
%         ylim([-0.01, 0.02]);
         xlim(100*[low_g,up_g])
%          xlim([-0.08,0.02])
       title('Nominal Interest Rate', 'fontsize', 15);
        %xlabel('g_t (Real Interest Rate Shock)', 'Fontsize', 12);
      xlabel('g_t (Aggregate Demand Shock)', 'Fontsize', 12);
      
       switch model
    case 'Taylor_pf' %|'Taylor_ZLB'
         legend('Active Regime', 'Passive Regime')  
       case 'Taylor_stc'
         legend('Active Regime', 'Passive Regime')        
       case 'Opt_pf' 
         legend('ZLB of Optimal in parfect forsight', 'no ZLB')
       case 'Opt_stc' 
         legend('ZLB of Optimal in stochastic', 'no ZLB')
      end  
       hold off
