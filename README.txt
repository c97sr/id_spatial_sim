Spatialsim source code
**********************

To do :
1 - Figure out why changing the max number of generations and the max number of infections
    in the parameter file doesn't seem to have an effect, even though it works for the 
    other parameters
2 - Get rid of all warnings ../src/libcpps/SR_Parameter.cpp:211:38: warning: non-constant array new length must be specified without parentheses around the type-id [-Wvla]
  paramlabels = new (string(*[dbgtmp]));
3 - Switch the random number generators from numerical recipes to gsl
