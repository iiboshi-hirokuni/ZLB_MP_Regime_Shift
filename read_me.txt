%%
%%    By Hirokuni Iiboshi (2016)
%% Å@"Monetary Policy Regime Shifts under the Zero Lower Bound:
%%    An application of a stochastic rational expectations equilibrium to a Markov switching DSGE model" 
%%    Economic Modelling, Vol. 52, p186-205  
%%
%%


Main programs are "Main_MS_Taylor_rule.m" and "Main_NK_Taylor_rule.m", the former of which calculates policy function of Markov switching version under ZLB, while the latter calculates those of a standard version under ZLB.

Values of parameters and options of ZLB and rational expectations types are selected from scripts of those matlab codes. And also, it is chosen that ranges of structral shocks.

These  options are located between line  15 and 50. Please set the values which you want.

The other options are numbers of grids such as shocks and forecasts by Gaussian Hermite, which are written in line 60 to 61.  And convergence accuracy is based on "Tolerance = 0.0025" in line 78. You can change those numbers.

When you excute these programs, then policy functions are saved in directories "output" for using making graphes.


For making graphes, there are 10 codes, in which each is named as "Plot_Fig_??" corresponding to Fig number of my paper. If you excute each of them, then a graph is made using  policy functions saved in directories "output".

Enjoy.
  