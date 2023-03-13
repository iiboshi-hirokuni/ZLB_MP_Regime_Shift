function [x,w] = GaussHermite(mu, sigma2, N)
  
HPoly1 = [ 1/pi^0.25 ];
HPoly2 = [ sqrt(2)/pi^0.25, 0];

for j = 1:N-1
    HPoly3 = [ sqrt(2/(j+1))*HPoly2, 0] - [0, 0, sqrt(j/(j+1))*HPoly1];
    HPoly1 = HPoly2;
    HPoly2 = HPoly3;
end

x1 = roots(HPoly3);

w1 = zeros(N,1);

for i = 1:N
    w1(i) = 1/(N)/(polyval(HPoly1,x1(i)))^2;
end

[ x, index] = sort(x1*sqrt(2*sigma2)+mu);
 
 w = w1(index)/sqrt(pi);