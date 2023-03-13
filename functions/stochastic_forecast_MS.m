% function  pi_1 = stochastic_forecast(PF) 

% global u g u_1 g_1 sigma_u sigma_g M

  [u_x,u_w] = GaussHermite(u_f, sigma_u^2, M);
  [g_x,g_w] = GaussHermite(g_f, sigma_g^2, M);
  
  coef = 1; % 0.5;
  
  if (min(u_x)<coef*mlow_u)|(min(g_x)<coef*mlow_g)|(max(u_x)>coef*mup_u)|(max(g_x)>coef*mup_g)
               pi_f0 = interp2(u,g,old_PF0_pi,u_f,g_f, 'cubic')+ pi_ss0;
               pi_f1 = interp2(u,g,old_PF1_pi,u_f,g_f, 'cubic')+ pi_ss1;                       
               
               y_f0 = interp2(u,g,old_PF0_y,u_f,g_f, 'cubic')+ y_ss0;
               y_f1 = interp2(u,g,old_PF1_y,u_f,g_f, 'cubic')+ y_ss1;               
                 
            % regime 0
                pi0_f = P(1,1)*pi_f0 +  P(2,1)*pi_f1;   
                y0_f = P(1,1)*y_f0 +  P(2,1)*y_f1; 
            % regime 1
                pi1_f = P(1,2)*pi_f0 +  P(2,2)*pi_f1;   
                y1_f = P(1,2)*y_f0 +  P(2,2)*y_f1;      
  else
    weight = kron(u_w, g_w);   
    pi0_ff = zeros(M,M);    
    y0_ff= zeros(M,M); 
    pi1_ff = zeros(M,M);    
    y1_ff= zeros(M,M); 
    
    
   for l = 1:M
       for m = 1:M
               pi_f0 = interp2(u,g,old_PF0_pi,u_x(l),g_x(m), 'cubic')+ pi_ss0;
               pi_f1 = interp2(u,g,old_PF1_pi,u_x(l),g_x(m), 'cubic')+ pi_ss1;                       
               
               y_f0 = interp2(u,g,old_PF0_y,u_x(l),g_x(m), 'cubic')+ y_ss0;
               y_f1 = interp2(u,g,old_PF1_y,u_x(l),g_x(m), 'cubic')+ y_ss1;               
                 
            % regime 0
               pi0_ff(l,m) = P(1,1)*pi_f0 +  P(2,1)*pi_f1;   
               y0_ff(l,m) = P(1,1)*y_f0 +  P(2,1)*y_f1; 
            % regime 1
               pi1_ff(l,m) = P(1,2)*pi_f0 +  P(2,2)*pi_f1;   
               y1_ff(l,m) = P(1,2)*y_f0 +  P(2,2)*y_f1;                   
      
%              pi_ff(l,m) = interp2(u,g,old_PF_pi,u_x(l),g_x(m), 'cubic');
%              y_ff(l,m) = interp2(u,g,old_PF_y,u_x(l),g_x(m), 'cubic');
       end
   end
  
  % ä˙ë“ÉCÉìÉtÉåÇÃçÏê¨ 
   pi0_f = sum( reshape(pi0_ff,M^2,1).*weight );
   pi1_f = sum( reshape(pi1_ff,M^2,1).*weight );
   
  % ä˙ë“_GDPÇÃçÏê¨ 
   y0_f = sum( reshape(y0_ff,M^2,1).*weight );
   y1_f = sum( reshape(y1_ff,M^2,1).*weight );
   
%   % perfect forsight
%   pi_f_1 = interp2(u,g,old_PF_pi,u_f,g_f, 'linear');
%   y_f_1 = interp2(u,g,old_PF_y,u_f,g_f, 'linear');
%     if (pi_f ~= pi_f_1) | (y_f ~= y_f_1)
%       [ pi_f  pi_f_1  y_f y_f_1 ]
%   end    
   
  end