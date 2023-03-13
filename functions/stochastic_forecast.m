% function  pi_1 = stochastic_forecast(PF) 

% global u g u_1 g_1 sigma_u sigma_g M

  [u_x,u_w] = GaussHermite(u_f, sigma_u^2, M);
  [g_x,g_w] = GaussHermite(g_f, sigma_g^2, M);
  
  coef = 1; % 0.5;
  
  if (min(u_x)<coef*mlow_u)|(min(g_x)<coef*mlow_g)|(max(u_x)>coef*mup_u)|(max(g_x)>coef*mup_g)
       pi_f = interp2(u,g,old_PF_pi,u_f,g_f, 'linear') + pi_ss; 
       y_f = interp2(u,g,old_PF_y,u_f,g_f, 'linear') + y_ss;      
  else
    weight = kron(u_w, g_w);   
    pi_ff = zeros(M,M); 
    y_ff= zeros(M,M); 
   for l = 1:M
       for m = 1:M
%              pi_ff(l,m) = interp2(u,g,old_PF_pi,u_x(l),g_x(m), 'linear');
%              y_ff(l,m) = interp2(u,g,old_PF_y,u_x(l),g_x(m), 'linear');
             
             pi_ff(l,m) = interp2(u,g,old_PF_pi,u_x(l),g_x(m), 'cubic') + pi_ss;
             y_ff(l,m) = interp2(u,g,old_PF_y,u_x(l),g_x(m), 'cubic') + y_ss;
       end
   end
  
  % Šú‘ÒƒCƒ“ƒtƒŒ‚Ìì¬ 
   pi_f = sum( reshape(pi_ff,M^2,1).*weight );
    
  % Šú‘Ò_GDP‚Ìì¬ 
   y_f = sum( reshape(y_ff,M^2,1).*weight );
  
%   % perfect forsight
%   pi_f_1 = interp2(u,g,old_PF_pi,u_f,g_f, 'linear');
%   y_f_1 = interp2(u,g,old_PF_y,u_f,g_f, 'linear');
%     if (pi_f ~= pi_f_1) | (y_f ~= y_f_1)
%       [ pi_f  pi_f_1  y_f y_f_1 ]
%   end    
   
  end