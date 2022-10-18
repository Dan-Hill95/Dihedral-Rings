% Copyright 2016 by Daniele Avitabile
% See https://www.maths.nottingham.ac.uk/personal/pmzda/
%
% If you use this code, please cite
% Daniele Avitabile, "Numerical computation of coherent structures in
% spatially-extended neural networks", Second International Conference on
% Mathematical Neuroscience, Antibes Juan-les-Pins, 2016

function branch = SecantContinuation(problem,u0,p0,stepperPars,name)

   %% Rename parameters
   s0            = stepperPars.s0;
   sMin          = stepperPars.sMin;
   sMax          = stepperPars.sMax;
   pMin          = stepperPars.pMin;
   pMax          = stepperPars.pMax;
   maxSteps      = stepperPars.maxSteps;
   nPrint        = stepperPars.nPrint;
   nSaveSol      = stepperPars.nSaveSol;
   iContPar      = stepperPars.iContPar;
   epsi          = stepperPars.finDiffEps;
   fsolveOptions = stepperPars.fsolveOptions;
   optNonlinIter = stepperPars.optNonlinIter;

      dataFolder    = stepperPars.dataFolder;
   if dataFolder(end)~='/'
     dataFolder = [dataFolder '/'];
   end

   PlotSolution    = stepperPars.PlotSolution;
   BranchVariables = stepperPars.BranchVariables;
   if isempty(BranchVariables) 
     numBranchVariables = 0;
   else
     numBranchVariables = size(BranchVariables(0,u0,p0),2);
   end
   if isempty(stepperPars.PlotBranchVariableId)
     iPlotBranchVariable = 4;
   else
     iPlotBranchVariable = 4 + stepperPars.PlotBranchVariableId;
   end

   ComputeEigenvalues = stepperPars.ComputeEigenvalues;
   numUnstableEigenvalues = -1;

   PlotSpectrum = stepperPars.PlotSpectrum;

   nDim          = size(u0,1);
   iU            = 1:nDim;
   iP            = nDim+1;

   %% Converge initial guess
   disp('*********** CONVERGE INITIAL GUESS *************');
   fsolveOptionsConverge = optimset(fsolveOptions,'Display','iter');
   [u0,fval,exitflag,output] = fsolve( @(u) problem(u,p0), u0, fsolveOptionsConverge );
   if exitflag <=0
     error('Failed to converge initial guess');
   end

   %% Perform a poorman step to launch secant continuation
   disp('*********** CONVERGE POOR MAN STEP *************');
   p1 = p0; p1(iContPar) = p0(iContPar) + s0/(10*sqrt(nDim));
   [u1,fval,exitflag,output] = fsolve( @(u) problem(u,p1), u0, fsolveOptionsConverge );
   if exitflag <=0
     error('Failed to converge poorman continuation initial step');
   end

   %% Initialise continuation
   step = 0; s = abs(s0);
   continuationFailed = false;
   reachedBound       = false;
   v0   = [u0; p0(iContPar)]; 
   v1   = [u1; p1(iContPar)];
   p    = p0;
   if ~isempty(ComputeEigenvalues)
     [W,D] = ComputeEigenvalues(v0(iU),p0);
     d = diag(D);
     d=d/abs(p(1));
     numUnstableEigenvalues = numel( find( real(d) > 0 ) );
   end
   DisplayStep(step,u0,p0,'init');
   SaveSolution(step,u0,p0,'init');
   SaveBranch(step,u0,p0,'init');
   bifDiagFigure = PlotBranch(branch,iPlotBranchVariable,[]);
   if ~isempty(PlotSolution)
     solPlotFigure = PlotSolution(u0,p0,[]);
   end
   if ~isempty(PlotSpectrum)
     solPlotSpectrum = PlotSpectrum(d,p0,[]);
   end

   %% Main continuation loop
   while step < maxSteps && ~continuationFailed && ~reachedBound

     % Predictor
     secant = (v1 - v0)/norm(v1 - v0, 2);
     v = v1 + secant * s;

     % Corrector
%      [v,fval,exitflag,output]  = fsolve( @(v) SecantCorrector(problem,v), v, fsolveOptions );
        Step = 1;
        while Step < 30
        u0=v;
        [F0,J0]=SecantCorrector(problem,u0);
        w0=-J0\F0;
        v=u0+w0;
        exitflag = 1;
        for k=1:p0(5)
        if norm(v(1+(k-1)*1000:k*1000)-u0(1+(k-1)*1000:k*1000))>1e-9 
            exitflag=0;
        end
        end
        if exitflag >0
        output.iterations=Step;
        break;
        else
        Step=Step+1;
        end
        end
     % Successful step
     if exitflag > 0 
   
       % Step control
       xi = optNonlinIter/output.iterations; 
       if xi < 0.5
	 xi = 0.5;
       elseif xi > 2;
	 xi = 2;
       end
       s = xi*s;
       %s = xi*s/sqrt(nDim);

       % Step limited by sMin and sMax
       if s > sMax
	 s = sMax;
       elseif s < sMin
	 s = sMin;
       end

       % Book keeping
       p(iContPar) = v(iP);
       step = step + 1;
       v0 = v1; 
       v1 = v;

       % Check if we hit the boundary
       if p(iContPar) < pMin || p(iContPar) > pMax
	 reachedBound = true;
       end

       %% Eventually compute eigenvalues
       if ~isempty(ComputeEigenvalues)
	 [W,D] = ComputeEigenvalues(v(iU),p);
	 d = diag(D);
     d=d/abs(p(1));
	 numUnstableEigenvalues = numel( find( real(d) > 0 ) );
       end

       % Output
       if mod(step,nPrint) == 0
	 DisplayStep(step,v(iU),p,[]);
	 PlotBranch(branch,iPlotBranchVariable,bifDiagFigure);
	 SaveBranch(step,v(iU),p,[]);
	 if ~isempty(PlotSolution)
	   PlotSolution(v(iU),p,solPlotFigure);
	 end
	 if ~isempty(PlotSpectrum)
	   solPlotSpectrum = PlotSpectrum(d,p,solPlotSpectrum);
	 end
       end
       if mod(step,nSaveSol) == 0
	 SaveSolution(step,v(iU),p,[]);
       end

     % Unsusccessful step
     else

       % Halving continuation step
       s = 0.5*s;
       disp(sprintf('Halving continuation step, s=%10.5e',s));

       % Stop if step is too small
       if s < sMin
	 continuationFailed = true;
       end

     end

     %% Output message
     if ~(step < maxSteps) 
       disp('Continuation ended: reached maximum number of steps');
     elseif reachedBound
       disp('Continuation ended: reached bound');
     elseif continuationFailed 
       disp(sprintf('Continuation failed: reached minimum step, s=%10.5e',s));
     end

   end

   function [F,J] = SecantCorrector(problem,v)

     % Set parameter
     p(iContPar) = v(iP);

     % Allocate
     F = zeros(size(v)); 
     if nargout > 1
       J = sparse(nDim,nDim); 
     end

     % Right-hand side
     if nargout > 1
       [F(iU),J(iU,iU)] = problem(v(iU),p);
     else
       F(iU) = problem(v(iU),p);
     end
     F(iP) = secant' * (v - v1) - s;

     % Jacobian
     if nargout > 1
       p(iContPar) = v(iP) + epsi;
       J(:,iP) = ( problem(v(iU),p) - F(iU) ) / epsi;
       J(iP,:) = secant';
       J = sparse(J);
     end

   end

   function DisplayStep(step,u,p,status)

     if strcmp(status,'init') 
       outString = '\n *********** START SECANT CONTINUATION *************\n';
       outString = [outString sprintf('%9s %5s %16s %14.6s %14.4s','STEP','STAB','PAR','2-NORM')];
       for i = 1:numBranchVariables
	 outString = [outString sprintf('%14.4s',['V(' num2str(i) ')'])];
       end
       outString = [outString '\n'];
       fprintf(outString);
     end

     formatString = '%9d %5d %16d %14.4e';
     for i = 1:numBranchVariables 
       formatString = [formatString ' %14.4e'];
     end
     formatString = [formatString '\n'];

     if numBranchVariables == 0
       outString = sprintf(formatString, step, numUnstableEigenvalues, p(iContPar), norm(u));
     else
       outString = sprintf(formatString, step, numUnstableEigenvalues, p(iContPar), norm(u), BranchVariables(step,u,p));
     end
     fprintf(outString);

   end

   function SaveSolution(step,u,p,status)

     if strcmp(status,'init') 
       if exist(dataFolder,'dir')
	 rmdir(dataFolder,'s');
       end
       mkdir(dataFolder);
     end

     fileName = [dataFolder sprintf('solution_%07d.mat', step)];
     if isempty(ComputeEigenvalues)
       save(fileName,'u','p');
     else
       save(fileName,'u','p','d','W');
     end

   end


   function SaveBranch(step,u,p,status)

     extraVariables = ~isempty(BranchVariables);

     if strcmp(status,'init') 
       if numBranchVariables == 0
	 branch = [step numUnstableEigenvalues p(iContPar) norm(u)];
       else
	 branch = [step numUnstableEigenvalues p(iContPar) norm(u) BranchVariables(step,u,p)];
       end
     else
       if numBranchVariables == 0
	 branch = [branch; step numUnstableEigenvalues p(iContPar) norm(u)];
       else
	 branch = [branch; step numUnstableEigenvalues p(iContPar) norm(u) BranchVariables(step,u,p)];
       end
     end
     fileName = [dataFolder 'branch.mat'];
     save(fileName,'branch');

   end

   function plotHandle = PlotBranch(branch,idVar,parentHandle)

     if isempty(parentHandle)
       scrsz = get(0,'ScreenSize');
       plotHandle = figure('Position',[scrsz(3)/4 scrsz(4)/2 scrsz(3)/4 scrsz(4)/4]);
       parentHandle = plotHandle;
     end

     figure(parentHandle)
     cla, hold on;
     plot(branch(:,3),branch(:,idVar),'b-')
     plot(branch(end,3),branch(end,idVar),'rd','MarkerSize',10);
     if ~isempty(ComputeEigenvalues)
       iUnstab = find(branch(:,2) > 0);
       %iStab   = setdiff(1:size(branch,1), iUnstab);
       iStab   = find(branch(:,2) == 0);
       plot(branch(iStab,3),branch(iStab,idVar),'b.');
       plot(branch(iUnstab,3),branch(iUnstab,idVar),'r.');
     else
       plot(branch(:,3),branch(:,idVar),'k.')
     end
     drawnow;

     print ('-dtiff', name)

   end

end
