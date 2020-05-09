**R code:**

**Script: postProcessing.R**

Source to produce .csv of summary stats for a set of simulations in working directory.

**Script: networkProcessing.R**

**RPdistributions(Node_connection_within_household, Node_connection_between_household, fname)**

Produces .csv of commuting distributions

**RPHHnet(Node_connection_within_household, Node_connection_between_household, fname)**

Produces .csv edge lists for HH-only networks from Rdata files.

**RPcombine(Node_connection_within_household, Node_connection_between_household)**

Produces single adjacency matrix for whole network as dgCMatrix.

**RPmcluster(G, N)**

Produces clustering coefficients up to order N for every node in network represented by dgCMatrix G.

**RPmclusterSample(G, mmax, sampleNumber)**

Produces clustering coefficients up to order mmax for sampleNumber nodes in network represented by adjacency dgCMatrix G. Note: this algorithm corresponds to DHCMsampleSpeed in MATLAB/octave code. Can be modified to produce summary statistics.

**mNbr(G, vx1, mMore, vi)**

Required for RPmclusterSample.

**MATLAB/Octave code:**

**RPcombine(G)**

Input edge list (columns I,J). Produces sparse adjacency matrix for a simple graph (binary, symmetric, no self-loops).

**DHCMall(G,m)**

Calculates clustering coefficients up to order m for every node in network represented by G. Produces summary statistics as required by plotting functions.

**DHCMsampling(G,m,sampleNumber)**

Calculates clustering coefficients up to order m for sampleNumber nodes in network represented by G. Produces summary statistics as required by plotting functions. If using Octave, statistics package is required (“pkg load statistics”).

**DHCMsamplingSpeed(G,m,sampleNumber)**

Calculates clustering coefficients up to order m for sampleNumber nodes in network represented by G. Produces summary statistics as required by plotting functions. This algorithm is slower on small networks but much faster on larger ones. If using Octave, statistics package is required (“pkg load statistics”).

**RpdiffAlphaProcessNets**

Set to produce one network per alpha value from C code output (one simulation per alpha). Can be modified to input/loop through different parameter sets.

**RpdiffAlphaProcessEpis**

Set to produce f=peak size, g=final size, h=cell array of epidemic curves, one row per alpha and one column per run, from C code output for a given R* (Rs) and R0 range. X is a matrix (one row per alpha, one column per run) of R0 values obtained from R code “postProcessing.R”. Can be modified to input/loop through different parameter sets.

**RpmclusPlotsError(cellMean)**

Plots CCm against alpha where cellMean is a 1 by a cell (a=number of alpha values). Each cell entry is a set of Cm summary statistics. If plotting for more than one network, is **RPcellMean.m** to calculate mean values. Max clustering order is calculated from cellMean{1}. All cell entries must be same size (i.e. to same order). **This was used to produce figure 4A.**

**RpcellMean.m**

Calculates pointwise mean of corresponding entries in 3 cell arrays. Can be modified to calculate for arbitrary number of cells.

RpincGen2(yCell,cCell)

Plots peak size/final size against clustering of each order, with linear fit and 95% confidence intervals. Inputs are n by 1 (or 1 by n) cell arrays with one entry per network. yCell{i}=[p,z] where p is a column vector of mean peak sizes and z a column vector of mean final sizes. cCell{i} is a 1 by a cell (a=number of alpha values), and each cell entry is a set of Cm summary statistics. **This was used to produce figures 4B and 4C.**

**Note: plotting figures 1 and 2 is a straightforward task from the outputs of postProcessing.R**
