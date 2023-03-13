clear

x = -4:0.01:4;

y = 1/sqrt(2*pi)*exp(-x.^2/2);

mu=0;
sigma2 =1;
N = 27;
[e,w] = GaussHermite(mu, sigma2, N);


[e w]
y1= zeros(size(x,2),1);

for j = 1:size(x,2)
for i = 1:size(e,1)
  if i < size(e,1)/2  
      if (x(j)> e(i)-1/2*abs(e(i)-e(i+1))) &(x(j)< e(i)+1/2*abs(e(i)-e(i+1)))
          y1(j)= w(i)/(abs(e(i)-e(i+1)));
      end
  else
      if (x(j)> e(i)-1/2*abs(e(i)-e(i-1))) &(x(j)< e(i)+1/2*abs(e(i)-e(i-1)))
          y1(j)= w(i)/(abs(e(i)-e(i-1)));
      end  
  end  
end
end

hold on
 area(x,y1,'FaceColor',[0.2 1 1],'LineStyle','none' )
 plot(x,y, 'r-','LineWidth',2.5)
 legend( 'Gaussian-Hermite','Normal Density')
 hold off
