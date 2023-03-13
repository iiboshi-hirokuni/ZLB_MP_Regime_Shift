%%
%%    By Hirokuni Iiboshi (2016)
%% 　"Monetary Policy Regime Shifts under the Zero Lower Bound:
%%    An application of a stochastic rational expectations equilibrium to a Markov switching DSGE model" 
%%    Economic Modelling, Vol. 52, p186-205  
%%
%%

clear all

% setting regime 0 or 1
regime = 1;   % 
picture = 1;  % draw IRF -> 1,  no drawing ->0

if picture == 1
   nsim = 100; % # of simulation
else   
   nsim = 100000;   % # of simulation
end
    
shock_u = 1;
shock_g = 1;

% Mixed polcy 
 p_00 = 0.7; p_11 = 0.7;
%  p_00_2 = 0.2; p_11_2 = 0.8;


%%  setting  type of model
                 model = 'Taylor_stc'
%              model = 'Taylor_pf'


    sigma_u = (0.154/100)*1;
    sigma_g = (1.524/100)*2;  %1.524
% % 
  rho_u = 0 ;
  rho_g = 0.8;
% 

% setting ranges of structural shocks
mlow_u = -0.1; mup_u = 0.1;  % markup shock
N_u = round((mup_u - mlow_u)/0.005);

mlow_g = -0.4; mup_g = 0.4; % real rate shock
N_g = round((mup_g - mlow_g)/0.002);


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
     
     % Aggressive Regime
     PF0_y = PF0_y_Taylor_pf;
     PF0_i = PF0_i_Taylor_pf;
     PF0_pi = PF0_pi_Taylor_pf;
     u = Taylor_pf_u;
     g = Taylor_pf_g;
   
     % Passive Regime
     PF1_y = PF1_y_Taylor_pf;
     PF1_i = PF1_i_Taylor_pf;
     PF1_pi = PF1_pi_Taylor_pf;
     u1 = Taylor_pf_u;
     g1 = Taylor_pf_g;
     
     load('./output/Taylor_pf_ZLB_no_switch.mat', 'PF_y_Taylor_pf');
     load('./output/Taylor_pf_ZLB_no_switch.mat', 'PF_i_Taylor_pf');
     load('./output/Taylor_pf_ZLB_no_switch.mat', 'PF_pi_Taylor_pf');         
     load('./output/Taylor_pf_ZLB_no_switch.mat', 'Taylor_pf_u');
     load('./output/Taylor_pf_ZLB_no_switch.mat', 'Taylor_pf_g');

   
     PF_y = PF_y_Taylor_pf;
     PF_i = PF_i_Taylor_pf;
     PF_pi = PF_pi_Taylor_pf;
%      u = Taylor_pf_u;
%      g = Taylor_pf_g;
%      
%  
     
 
 case 'Taylor_stc' 
     load('./output/Taylor_stc_ZLB.mat', 'PF0_y_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF0_i_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF0_pi_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF1_y_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF1_i_Taylor_stc');
     load('./output/Taylor_stc_ZLB.mat', 'PF1_pi_Taylor_stc');     
     load('./output/Taylor_stc_ZLB.mat', 'Taylor_stc_u');
     load('./output/Taylor_stc_ZLB.mat', 'Taylor_stc_g');

%      if regime == 1
     % Aggressive Regime
     PF0_y = PF0_y_Taylor_stc;
     PF0_i = PF0_i_Taylor_stc;
     PF0_pi = PF0_pi_Taylor_stc;
     u = Taylor_stc_u;
     g = Taylor_stc_g;
%      end
     
     % Passive Regime 
     PF1_y = PF1_y_Taylor_stc;
     PF1_i = PF1_i_Taylor_stc;
     PF1_pi = PF1_pi_Taylor_stc;
     u1 = Taylor_stc_u;
     g1 = Taylor_stc_g;
   
%   case 'Taylor_stc_no_switch' 
 
      load('./output/Taylor_stc_ZLB_no_switch.mat', 'PF_y_Taylor_stc');
     load('./output/Taylor_stc_ZLB_no_switch.mat', 'PF_i_Taylor_stc');
     load('./output/Taylor_stc_ZLB_no_switch.mat', 'PF_pi_Taylor_stc');         
     load('./output/Taylor_stc_ZLB_no_switch.mat', 'Taylor_stc_u');
     load('./output/Taylor_stc_ZLB_no_switch.mat', 'Taylor_stc_g');

   
     PF_y = PF_y_Taylor_stc;
     PF_i = PF_i_Taylor_stc;
     PF_pi = PF_pi_Taylor_stc;
%      u = Taylor_stc_u;
%      g = Taylor_stc_g;
     
%      case 'Taylor_pf_no_switch' 
%  

%      
 end    
        

u_sim = 0;;  % マークアップショック
g_sim = 0;      % 実質金利ショック
pi_sim = 0; 
y_sim = 0; %  
i_sim = 0; %
st_save = 0;

sub_y = -30;

pi0 = interp2(u, g, PF0_pi, 0, 0, 'linear');
y0  = interp2(u, g, PF0_y, 0, 0, 'linear');
i0  = interp2(u, g, PF0_i, 0, 0, 'linear');

st=0;

sd_u = 50;
rand('seed',sd_u);
stt = rand(1,nsim);

sd = 1000;
randn('seed',sd);
shock = randn(2,nsim);

for j = 2:nsim  
  d=0;  
    
      u_sim(j) = 0;
      g_sim(j) = 0;
      
      if shock_u == 1
        u_sim(j) = rho_u * u_sim(j-1)+sigma_u*shock(1,j);
      end
      if shock_g == 1
        g_sim(j) = rho_g * g_sim(j-1)+sigma_g*shock(2,j);
      end
      
      if st==0;
          p = 1-p_00;  % prob(st=1 | st=0)
      else
          p = p_11;    % prob(st=1 | st=1)
      end
      
      %       st = binornd(1,p);
      
      if stt(j) < p
       st = 1;
      else
       st = 0;
      end 
      
      st_save(j) = st;
      
  if st == 0;           
      pi_sim(j) = (interp2(u, g, PF0_pi, u_sim(j), g_sim(j), 'linear')-pi0)*100;
      y_sim(j) = interp2(u, g, PF0_y, u_sim(j), g_sim(j), 'linear')*100;
      i_sim(j) = interp2(u, g, PF0_i, u_sim(j), g_sim(j), 'linear')*100;
  else
      pi_sim(j) = (interp2(u, g, PF1_pi, u_sim(j), g_sim(j), 'linear')-pi0)*100;
      y_sim(j) = interp2(u, g, PF1_y, u_sim(j), g_sim(j), 'linear')*100;
      i_sim(j) = interp2(u, g, PF1_i, u_sim(j), g_sim(j), 'linear')*100;
  end      
      
   while d==0  
   if isnan(pi_sim(j))==1|isnan(i_sim(j))==1|isnan(y_sim(j))==1 | y_sim(j)<sub_y
        u_sim(j) = rho_u * u_sim(j-1)+sigma_u*randn(1,1);
       g_sim(j) = rho_g * g_sim(j-1)+sigma_g*randn(1,1);
       
     if st == 0;        
       pi_sim(j) = (interp2(u, g, PF0_pi, u_sim(j), g_sim(j), 'linear')-pi0)*100;
       y_sim(j) = interp2(u, g, PF0_y, u_sim(j), g_sim(j), 'linear')*100;
       i_sim(j) = interp2(u, g, PF0_i, u_sim(j), g_sim(j), 'linear')*100;
     else        
        pi_sim(j) = (interp2(u, g, PF1_pi, u_sim(j), g_sim(j), 'linear')-pi0)*100;
        y_sim(j) = interp2(u, g, PF1_y, u_sim(j), g_sim(j), 'linear')*100;
        i_sim(j) = interp2(u, g, PF1_i, u_sim(j), g_sim(j), 'linear')*100;
     end  
  else
      d=1;  
   end 
   end
  
end

mean_y = num2str(mean(y_sim)); std_y = num2str(std(y_sim));
mean_pi = num2str(mean(pi_sim)); std_pi = num2str(std(pi_sim));
mean_i = num2str(mean(i_sim)); std_i = num2str(std(i_sim));

switch model
        case 'Taylor_stc';
%            display(['']);
             display('Mixed Policy');
        case 'Taylor_pf';
           display(['']);
             display('aggressive regime');
        case 'Taylor_pf_no_switch'   
            display(['']);
            display('Fixed Regime');
        case 'Taylor_stc_no_switch'   
%             display(['']);
            display('\n Fixed Regime');
end
display(['nsim=' num2str(nsim) ]);
display(['mean_y=' mean_y ', std_y=' std_y ])
display(['mean_pi=' mean_pi ', std_pi=' std_pi ])
display(['mean_i=' mean_i ', std_i=' std_i ])

if  picture == 1
figure(1)
subplot(3,3,1); hist(y_sim,15); title('y','FontSize',14 ); xlim([-15 15]);
subplot(3,3,5); hist(pi_sim,15); title('\pi','FontSize',14 ); xlim([-1.5 1]);
subplot(3,3,9); hist(i_sim,15); title('i','FontSize',14 ); xlim([-1.5 1.5]);
subplot(3,3,4); scatter(y_sim,pi_sim,'.'); title('\pi vs y','FontSize',14 ); xlim([-15 15]); ylim([-1.5 1]);
subplot(3,3,7); scatter(y_sim,i_sim,'.'); title('i vs y','FontSize',14 ); xlim([-15 15]);
subplot(3,3,8); scatter(pi_sim,i_sim,'.'); title('i vs \pi','FontSize',14 ); xlim([-1.5 1]);
subplot(3,3,2); scatter(pi_sim,y_sim,'.'); title('y vs \pi','FontSize',14 ); xlim([-1.5 1]); ylim([-15 15]);
subplot(3,3,3); scatter(i_sim,y_sim,'.'); title('y vs i','FontSize',14 ); xlim([-1.5 1.5]); ylim([-15 15]);
subplot(3,3,6); scatter(i_sim,pi_sim,'.'); title('pi vs i','FontSize',14 ); xlim([-1.5 1.5]); ylim([-1.5 1]);
end

u1_sim = 0;;     % markup shock
g1_sim = 0;      % real rate shock
pi1_sim = 0; 
y1_sim = 0; %  
i1_sim = 0; %

pi1 = interp2(u1, g1, PF1_pi, 0, 0, 'linear');
st = 0;

for j = 2:nsim
  
  d=0;  
    
     u1_sim(j) = 0;
     g1_sim(j) = 0;

     if shock_u == 1
       u1_sim(j) = rho_u * u1_sim(j-1)+sigma_u*shock(1,j);
     end
     if shock_g == 1
       g1_sim(j) = rho_g * g1_sim(j-1)+sigma_g*shock(2,j);
     end
     
%       if st==0;
%           p = 1-p_00_2;  % prob(st=1 | st=0)
%       else
%           p = p_11_2;    % prob(st=1 | st=1)
%       end
%       
%       st = binornd(1,p);
%       
%   if st == 0;           
      pi1_sim(j) = (interp2(u1, g1, PF_pi, u1_sim(j), g1_sim(j), 'linear')-pi1)*100;
      y1_sim(j) = interp2(u1, g1, PF_y, u1_sim(j), g1_sim(j), 'linear')*100;
      i1_sim(j) = interp2(u1, g1, PF_i, u1_sim(j), g1_sim(j), 'linear')*100;
%   else
%       pi1_sim(j) = (interp2(u, g, PF1_pi, u1_sim(j), g1_sim(j), 'linear')-pi0)*100;
%       y1_sim(j) = interp2(u, g, PF1_y, u1_sim(j), g1_sim(j), 'linear')*100;
%       i1_sim(j) = interp2(u, g, PF1_i, u1_sim(j), g1_sim(j), 'linear')*100;
%   end      
       
  
    
   while d==0  
    if isnan(pi1_sim(j))==1|isnan(i1_sim(j))==1|isnan(y1_sim(j))==1 | y1_sim(j)<sub_y
       u1_sim(j) = rho_u * u1_sim(j-1)+sigma_u*randn(1,1);
       g1_sim(j) = rho_g * g1_sim(j-1)+sigma_g*randn(1,1);
       
%        if st == 0;        
       pi1_sim(j) = (interp2(u, g, PF_pi, u1_sim(j), g1_sim(j), 'linear')-pi0)*100;
       y1_sim(j) = interp2(u, g, PF_y, u1_sim(j), g1_sim(j), 'linear')*100;
       i1_sim(j) = interp2(u, g, PF_i, u1_sim(j), g1_sim(j), 'linear')*100;
%      else        
%         pi1_sim(j) = (interp2(u, g, PF1_pi, u1_sim(j), g1_sim(j), 'linear')-pi0)*100;
%         y1_sim(j) = interp2(u, g, PF1_y, u1_sim(j), g1_sim(j), 'linear')*100;
%         i1_sim(j) = interp2(u, g, PF1_i, u1_sim(j), g1_sim(j), 'linear')*100;
%      end  
   else
      d=1;  
   end 
   end
  
end

mean_y1 = num2str(mean(y1_sim)); std_y1 = num2str(std(y1_sim));
mean_pi1 = num2str(mean(pi1_sim)); std_pi1 = num2str(std(pi1_sim));
mean_i1 = num2str(mean(i1_sim)); std_i1 = num2str(std(i1_sim));

% display('');
display('Fixed Policy');
display(['mean_y1=' mean_y1 ', std_y1=' std_y1 ])
display(['mean_pi1=' mean_pi1 ', std_pi1=' std_pi1 ])
display(['mean_i1=' mean_i1 ', std_i1=' std_i1 ])

if  picture == 1
figure('Name','ZLB')  

subplot(3,3,1); 
[density,x1]=ksdensity(y_sim);
[density,x2]=ksdensity(y1_sim);
   hold on
        plot(x1,density,'LineStyle','-','Color','b',...
        'LineWidth',2.5);  
       plot(x2,density,'LineStyle',':','Color','r',...
        'LineWidth',2.5); title('y','FontSize',14 ) ;
    xlim([-10 10]);
%     switch model
%         case 'Taylor_stc';
%            legend('Aggressive Regime','Passive Regime' )
%         case 'Taylor_pf';
%            legend('Aggressive Regime','Passive Regime' )
%         case 'Taylor_pf_no_switch'   
%             legend('Fixed Regime','Passive Regime')
%         case 'Taylor_stc_no_switch'   
%             legend('Fixed Regime','Passive Regime' )
%     end
   legend('Mixed Policy','Fixed Policy' )
   
   hold off 
subplot(3,3,5); [density,x1]=ksdensity(pi_sim);
        [density,x2]=ksdensity(pi1_sim);
   hold on
        plot(x1,density,'LineStyle','-','Color','b',...
        'LineWidth',2.5);  
       plot(x2,density,'LineStyle',':','Color','r',...
        'LineWidth',2.5); title('\pi','FontSize',14 ) ; xlim([-0.5 0.5]);
   hold off 
subplot(3,3,9); [density,x1]=ksdensity(i_sim);
       [density,x2]=ksdensity(i1_sim);
   hold on
        plot(x1,density,'LineStyle','-','Color','b',...
        'LineWidth',2.5);  
       plot(x2,density,'LineStyle',':','Color','r',...
        'LineWidth',2.5); title('i','FontSize',14 ) ;
     xlim([-2 2]); xlabel('i','FontSize',14);
   hold off 
subplot(3,3,4); hold on
   scatter(y_sim,pi_sim,'b.');
   scatter(y1_sim,pi1_sim,'r.');
%    scatter(y_sim,pi_sim,'b.');
   
   ylabel('\pi','FontSize',14);
    title('\pi vs y','FontSize',14 ); xlim([-10 10]); ylim([-0.5 0.5]); hold off
subplot(3,3,7); hold on
   scatter(y_sim,i_sim,'b.');  
   scatter(y1_sim,i1_sim,'r.');
%    scatter(y_sim,i_sim,'b.');  
   
   xlabel('y','FontSize',14);  ylabel('i','FontSize',14);
   title('i vs y','FontSize',14 ); ylim([-2 2]); xlim([-10 10]);  hold off
subplot(3,3,8); hold on
   scatter(pi_sim,i_sim,'b.'); 
  scatter(pi1_sim,i1_sim,'r.');
%    scatter(pi_sim,i_sim,'b.');   
   xlabel('\pi','FontSize',14);
   title('i vs \pi','FontSize',14 ); ylim([-2 2]);xlim([-0.5 0.5]); hold off
subplot(3,3,2); hold on
%   scatter(pi_sim,y_sim,'b.'); 
  scatter(pi1_sim,y1_sim,'r.');
  scatter(pi_sim,y_sim,'b.'); 
  
  title('y vs \pi','FontSize',14 );  
  xlim([-0.5 0.5]); ylim([-10 10]);  hold off
subplot(3,3,3); hold on
%     scatter(i_sim,y_sim,'b.');
  scatter(i1_sim,y1_sim,'r.');
  scatter(i_sim,y_sim,'b.'); 
  
  title('y vs i','FontSize',14 );    
  xlim([-2 2]); ylim([-10 10]);  hold off
subplot(3,3,6); hold on
%   scatter(i_sim,pi_sim,'b.'); 
  scatter(i1_sim,pi1_sim,'r.');
  scatter(i_sim,pi_sim,'b.'); 
  
  title('\pi vs i','FontSize',14 );  
  xlim([-2 2]); ylim([-0.5 0.5]);  hold off
  
  

figure
subplot(3,1,1)
 hold on
     hh(1)= area(st_save*1000 );
     hh(2)=  area(st_save*-1000);
       set(hh(1),'LineStyle','none')
       set(hh(1),'FaceColor',[0.7 1 1])
       set(hh(2),'LineStyle','none')
       set(hh(2),'FaceColor',[0.7 1 1])
       
     plot(y_sim,'b-','LineWidth',2.5);
     plot(y1_sim,'r:','LineWidth',2.5);
    
     legend('Passive Regime', '', 'Regime Swiching Policy','Fixed Regime Policy' )
     title('Output','FontSize',14 );
     ylim([-20 15]);
 hold off  
subplot(3,1,2)
  hold on
      hh(1)= area(st_save*1000 );
     hh(2)=  area(st_save*-1000);
       set(hh(1),'LineStyle','none')
       set(hh(1),'FaceColor',[0.7 1 1])
       set(hh(2),'LineStyle','none')
       set(hh(2),'FaceColor',[0.7 1 1])
       
     plot(pi_sim,'b-','LineWidth',2.5);
     plot(pi1_sim,'r:','LineWidth',2.5);
     legend('Passive Regime', '', 'Regime Swiching Policy','Fixed Regime Policy' )
     title('Inflation','FontSize',14 );
      ylim([-0.7 0.6]);
 hold off  
subplot(3,1,3)
    hold on
      hh(1)= area(st_save*1000 );
     hh(2)=  area(st_save*-1000);
       set(hh(1),'LineStyle','none')
       set(hh(1),'FaceColor',[0.7 1 1])
       set(hh(2),'LineStyle','none')
       set(hh(2),'FaceColor',[0.7 1 1])
     plot(i_sim,'b-','LineWidth',2.5);
     plot(i1_sim,'r:','LineWidth',2.5);
     legend('Passive Regime', '', 'Regime Swiching Policy','Fixed Regime Policy' )
     title('Nominal Interest Rate','FontSize',14 );
      ylim([-1.0 3]);
 hold off  
  
end 


% if  picture == 1
% figure(2)  
% subplot(3,3,1); 
% [density,x1]=ksdensity(y_sim);
% [density,x2]=ksdensity(y1_sim);
%    hold on
%         plot(x1,density,'LineStyle','-','Color','b',...
%         'LineWidth',2.5);  
%        plot(x2,density,'LineStyle',':','Color','r',...
%         'LineWidth',2.5); title('y','FontSize',14 ) ;
%     xlim([-10 5]);
%     
%     switch model
%         case 'Taylor_stc';
%            legend('Aggressive Regime','Passive Regime' )
%         case 'Taylor_pf';
%            legend('Aggressive Regime','Passive Regime' )
%         case 'Taylor_pf_no_switch'   
%             legend('Fixed Regime','Passive Regime' )
%         case 'Taylor_stc_no_switch'   
%             legend('Fixed Regime','Passive Regime' )
%     end
%    hold off 
% subplot(3,3,5); [density,x1]=ksdensity(pi_sim);
%         [density,x2]=ksdensity(pi1_sim);
%    hold on
%         plot(x1,density,'LineStyle','-','Color','b',...
%         'LineWidth',2.5);  
%        plot(x2,density,'LineStyle',':','Color','r',...
%         'LineWidth',2.5); title('pi','FontSize',14 ) ; xlim([-10 10]);
%    hold off 
% subplot(3,3,9); [density,x1]=ksdensity(i_sim);
%        [density,x2]=ksdensity(i1_sim);
%    hold on
%         plot(x1,density,'LineStyle','-','Color','b',...
%         'LineWidth',2.5);  
%        plot(x2,density,'LineStyle',':','Color','r',...
%         'LineWidth',2.5); title('i','FontSize',14 ) ;
%      xlim([-2 6]);
%    hold off 
% subplot(3,3,4); hold on
%    scatter(y1_sim,pi1_sim,'r.');
%    scatter(y_sim,pi_sim,'b.'); 
%     title('pi vs y','FontSize',14 ); xlim([-10 5]); ylim([-10 10]); hold off
% subplot(3,3,7); hold on
%    scatter(y1_sim,i1_sim,'r.');
%    scatter(y_sim,i_sim,'b.');   
%    title('i vs y','FontSize',14 ); ylim([-2 6]); xlim([-10 5]);  hold off
% subplot(3,3,8); hold on
%   scatter(pi1_sim,i1_sim,'r.');
%    scatter(pi_sim,i_sim,'b.');   
%    title('i vs pi','FontSize',14 ); ylim([-2 6]);xlim([-10 10]); hold off
% subplot(3,3,2); hold on
%   scatter(pi1_sim,y1_sim,'r.');
%   scatter(pi_sim,y_sim,'b.'); title('y vs pi','FontSize',14 );  
%   xlim([-10 10]); ylim([-10 5]);  hold off
% subplot(3,3,3); hold on
%   scatter(i1_sim,y1_sim,'r.');
%   scatter(i_sim,y_sim,'b.'); title('y vs i','FontSize',14 );  
%   xlim([-2 6]); ylim([-10 5]);  hold off
% subplot(3,3,6); hold on
%   scatter(i1_sim,pi1_sim,'r.');
%   scatter(i_sim,pi_sim,'b.'); title('pi vs i','FontSize',14 );  
%   xlim([-2 6]); ylim([-10 10]);  hold off
% end 
%  
% 
% pi0 = 100*interp2(u, g, PF_pi, 0, 0, 'linear')
% y0  = 100*interp2(u, g, PF_y, 0, 0, 'linear')
% i0  = 100*interp2(u, g, PF_i, 0, 0, 'linear')
% pi1 = 100*interp2(u1, g1, PF1_pi, 0, 0, 'linear')
% y1  = 100*interp2(u1, g1, PF1_y, 0, 0, 'linear')
% i1  = 100*interp2(u1, g1, PF1_i, 0, 0, 'linear')

