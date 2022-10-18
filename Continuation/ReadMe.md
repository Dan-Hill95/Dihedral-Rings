*These codes have either been written in 2021 by Dan J. Hill or adapted from codes by David J.B. Lloyd (University of Surrey) and Daniele Avitabile, "Numerical computation of coherent structures in spatially-extended neural networks", Second International Conference on Mathematical Neuroscience, Antibes Juan-les-Pins, 2016.*

**For a more thorough tutorial in secant continuation and access to the original codes, see https://github.com/danieleavitabile/continuation-spatially-extended-systems-tutorial.**

Rather than accessing any of the available folders, you should initialise this code by running

*>> Init.m*


**After running Init.m, the following guide should appear in your command window:**

These codes are for the continuation of localised dihedral patterns in the 2-3 Swift-Hohenberg equation.
Here are your next steps:

1. If you want to solve the (n+1)-dimensional algebraic matching condition for localised dihedral
   patterns, type:
 
*>> a = MatchSoln(x, m, r_max, mu);*
 
   where x=[x[0],x[1],...,x[n]] is your initial guess, m is the dihedral lattice,
   r_max is the radial boundary and mu gives an approximate localisation. For example:
 
*>> a = MatchSoln((1/10) * ones(1,11), 6, 100, 1e-4);*
 
2. If you want to find and continue localised dihedral patterns in the 2-3 SH equation, type:
 
*>> branch = Cont_Patch(x, p, Dir);*
 
   where p=[mu, gamma, kappa, m, n+1] contains the parameters of the system including the bifurcation
   parameter mu, the respective quadratic and cubic coefficients (gamma, kappa), and the dimension 
   (n+1) of the reduced ODE system. Dir indicates the direction and step size of the continuation
   routine, which must be 'pl' (plus), 'mn' (minus), 'sp' (small plus), or 'sm' (small minus).
   For example:
 
*>> branch = Cont_Patch([-1,2], [0.02, 1.6, -1, 2, 5], 'pl');*
 
3. If you want to plot the bifurcation curve of a localised solution, type:
 
*>> ExploreBifurcationDiagram('FolderName/branch.mat',idVar);*
 
   where FolderName is the data folder that is created by Cont_Patch, and idVar is the column of the
   branch measures to be plotted. For example, for a data folder called 'D2_Patch_pl', you would type:
 
*>> ExploreBifurcationDiagram('D2_Patch_pl/branch.mat',5);*
 
