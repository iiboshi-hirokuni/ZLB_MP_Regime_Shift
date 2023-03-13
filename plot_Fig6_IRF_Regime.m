%%
%%    By Hirokuni Iiboshi (2016)
%% Å@"Monetary Policy Regime Shifts under the Zero Lower Bound:
%%    An application of a stochastic rational expectations equilibrium to a Markov switching DSGE model" 
%%    Economic Modelling, Vol. 52, p186-205  
%%
%%

shock1 = 0.0;    % markup shock

shock_sign =0;  % postive=1, negative=0 
if shock_sign ==1
   shock2 = 0.1; %-0.1;    % real rate shocks
else
    shock2 = -0.1;
end    
   
rho_u = 0 ;  
rho_g = 0.8;

 % model = 'Taylor_pf'
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

     % Active regime
     PF0_y = PF0_y_Taylor_pf;
     PF0_i = PF0_i_Taylor_pf;
     PF0_pi = PF0_pi_Taylor_pf;
     u_ZLB = Taylor_pf_u;
     g_ZLB = Taylor_pf_g;
   
     % passive regime
     PF1_y = PF1_y_Taylor_pf;
     PF1_i = PF1_i_Taylor_pf;
     PF1_pi = PF1_pi_Taylor_pf;
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

     PF0_y = PF0_y_Taylor_stc;
     PF0_i = PF0_i_Taylor_stc;
     PF0_pi = PF0_pi_Taylor_stc;
     u = Taylor_stc_u;
     g = Taylor_stc_g;
  
     PF1_y = PF1_y_Taylor_stc;
     PF1_i = PF1_i_Taylor_stc;
     PF1_pi = PF1_pi_Taylor_stc;
     
     load('./output/Taylor_stc_ZLB_no_switch.mat', 'PF_y_Taylor_stc');
     load('./output/Taylor_stc_ZLB_no_switch.mat', 'PF_i_Taylor_stc');
     load('./output/Taylor_stc_ZLB_no_switch.mat', 'PF_pi_Taylor_stc'); 
     
     PF_y = PF_y_Taylor_stc;
     PF_i = PF_i_Taylor_stc;
     PF_pi = PF_pi_Taylor_stc;    
     
end


% u=u_ZLB; g=g_ZLB;

%% Impulse Response Functions
nsim = 100;    % # of Monte Calro
hh = 20;       % horizons of response  

%% Regime 0
u_sim1 = shock1; %-0.05; % markup shock
g_sim1 = shock2; %-0.1;    % real rate shock
pi_sim1 = 0; 
y_sim1 = 0; %  
i_sim1 = 0; %

 pi_0 = interp2(u, g, PF0_pi, 0, 0, 'linear');
  y_0 = interp2(u, g, PF0_y, 0, 0, 'linear');
  i_0 = interp2(u, g, PF0_i, 0, 0, 'linear');

for j = 1:nsim
  %   
  pi_sim1(j) = interp2(u, g, PF0_pi, u_sim1(j), g_sim1(j), 'linear')- pi_0;
  y_sim1(j) = interp2(u, g, PF0_y, u_sim1(j), g_sim1(j), 'linear')- y_0;
  i_sim1(j) = interp2(u, g, PF0_i, u_sim1(j), g_sim1(j), 'linear')- i_0;
  
  if j < nsim
    u_sim1(j+1) = rho_u * u_sim1(j);
    g_sim1(j+1) = rho_g * g_sim1(j);
  end
end


%% Regime  1

u_sim2 = u_sim1; %-0.05;  % markup shock
g_sim2 = g_sim1;      % real rate shock
pi_sim2 = 0; 
y_sim2 = 0; %  
i_sim2 = 0; %

  pi_1 = interp2(u, g, PF1_pi, 0, 0, 'linear');
  y_1 = interp2(u, g, PF1_y, 0, 0, 'linear');
  i_1 = interp2(u, g, PF1_i, 0, 0, 'linear');
  
for j = 1:nsim
  pi_sim2(j) = interp2(u, g, PF1_pi, u_sim2(j), g_sim2(j), 'linear')- pi_1  ;
  y_sim2(j) = interp2(u, g, PF1_y, u_sim2(j), g_sim2(j), 'linear')- y_1;
  i_sim2(j) = interp2(u, g, PF1_i, u_sim2(j), g_sim2(j), 'linear')- i_1;

  if j < nsim
    u_sim2(j+1) = rho_u * u_sim2(j);
    g_sim2(j+1) = rho_g * g_sim2(j);
  end
end

%% no Regime

u_sim3 = u_sim1; %-0.05;  % markup shock
g_sim3 = g_sim1;      % real rate shock
pi_sim3 = 0; 
y_sim3 = 0; %  
i_sim3 = 0; %

  pi_3 = interp2(u, g, PF_pi, 0, 0, 'linear');
  y_3 = interp2(u, g, PF_y, 0, 0, 'linear');
  i_3 = interp2(u, g, PF_i, 0, 0, 'linear');
  
for j = 1:nsim
  pi_sim3(j) = interp2(u, g, PF_pi, u_sim3(j), g_sim3(j), 'linear')- pi_3  ;
  y_sim3(j) = interp2(u, g, PF_y, u_sim3(j), g_sim3(j), 'linear')- y_3;
  i_sim3(j) = interp2(u, g, PF_i, u_sim3(j), g_sim3(j), 'linear')- i_3;

  if j < nsim
    u_sim3(j+1) = rho_u * u_sim3(j);
    g_sim3(j+1) = rho_g * g_sim3(j);
  end
end



figure
   y_sim = [0 y_sim1(1:hh);    0 y_sim2(1:hh); 0 y_sim3(1:hh) ]*100;
   pi_sim = [0 pi_sim1(1:hh); 0 pi_sim2(1:hh); 0 pi_sim3(1:hh) ]*100;
    i_sim = [0 i_sim1(1:hh);    0 i_sim2(1:hh); 0 i_sim3(1:hh) ]*100;
%    i_sim = [i_0 i_sim1(1:hh);    i_1 i_sim2(1:hh); i_3 i_sim3(1:hh) ]*100;
subplot(3,1,1)   
hold on
plot(0:hh, y_sim(1,:),'b-','LineWidth',2  );  % Active
 plot(0:hh, y_sim(2,:),'r--','LineWidth',2  ); % Passive
  plot(0:hh, y_sim(3,:),'c:','LineWidth',2  ); % Fixed
plot(0:hh, zeros(1,hh+1),'g:','LineWidth',2 );
% ylim([-15 2])
hold off
title('Output', 'fontsize', 15); 
legend('Active Regime', 'Passive Regime','Fixed Regime'); 
subplot(3,1,2)
hold on
plot(0:hh, pi_sim(1,:),'b-','LineWidth',2  );
 plot(0:hh, pi_sim(2,:),'r--','LineWidth',2  );
  plot(0:hh, pi_sim(3,:),'c:','LineWidth',2  );
plot(0:hh, zeros(1,hh+1),'g:','LineWidth',2 );
% ylim([-2 10])
hold off
title('Inflation', 'fontsize', 15); 
% legend('Active Regime', 'Passive Regime'); 
legend('Active Regime', 'Passive Regime','Fixed  Regime'); 
subplot(3,1,3)
hold on
plot(0:hh, i_sim(1,:),'b-','LineWidth',2  );
 plot(0:hh, i_sim(2,:),'r--','LineWidth',2  );
 plot(0:hh, i_sim(3,:),'c:','LineWidth',2  );
plot(0:hh, zeros(1,hh+1),'g:','LineWidth',2 );
% ylim([-0.2 6])
hold off
title('Nominal Interest Rate', 'fontsize', 15);
% legend('Active Regime', 'Passive Regime'); 
legend('Active Regime', 'Passive Regime','Fixed  Regime'); 