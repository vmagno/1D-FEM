1;

function [y, dy] = exactResult1(f, T, x)
  y = (f / (2*T)) .* x .* (1 - x);
  dy = (f / (2*T)) * (1 - 2 .*x);
endfunction

function [y, dy] = exactResult2(k, x)
  num = exp((x - 1) ./ k);
  denom = (exp(-1/k) - 1);
  y = (num - 1) ./ denom;
  dy = (1 / (k * denom)) * num;
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k, c, b and f are the parameters of the differential equation
%   ku'' + cu' + bu + f = 0
%

function errorNorm(nEls,nodeCoords,connect,elDof,dFreedom,pDeg,pType,u,k,c,b,f)
  %%%%%%%%%%%%%
  % Norm of error is Sum_elem ( int((uexact - uh)^2) )
  %
  % With gaussian quadrature: Sum_elem( Sum_gpoints ( weight * (uexact - uh)^2 * h/2 ) )

  %figure(2)
  x = 0:0.1:1;
  
  subplot(2,2,3)
  hold on
  axis square
  plot(x,exactResult2(k(x), x))
  
  for iElem = 1:nEls
    for i = 1:elDof(iElem)
        c = dFreedom(iElem,i);
        uL(iElem,i) = u(c);
    end
  end
  
  
  TotalError = 0;
  TotalDerivError = 0;
  [xiQ,wQ] = gQuad;

  for iElem = 1 : nEls
    x1 = nodeCoords(connect(iElem, 1));
    x2 = nodeCoords(connect(iElem, 2));
    h = x2 - x1;
    jac = 2 / h;
    
    nQ = 3; % order of gaussian quadrature
    for iG = 1:nQ;
    
      xi = xiQ(iG,nQ);
      wq = wQ(iG,nQ);
      xq = h/2*xi + (x1+x2)/2;   %global coordinate of xi
      
      % Get value of computed solution at point xi
      [N,dN] = shape(xi,iElem,pDeg,pType);
      uh = 0;
      duh = 0;
      for iDof = 1:elDof(iElem)
        uh = uh + N(iDof) * uL(iElem, iDof);
        duh = duh + dN(iDof) * uL(iElem, iDof) * 2 / h;
      end
      
      % Value of exact solution at point xi (xq)
      [uexact, duexact] = exactResult2(k(xq), xq);
      
      % Value at this gaussian point
      TotalError = TotalError + wq * (uexact - uh)^2 * h / 2;
      TotalDerivError = TotalDerivError + wq * (duexact - duh)^2 * h / 2;
      
    end
    
    
  end

  L2Norm = sqrt(TotalError)
  L2DerivNorm = sqrt(TotalDerivError)

endfunction