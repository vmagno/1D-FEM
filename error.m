1;

function [y, dy] = exactResult1(f, T, x)
  y = (f / (2*T)) .* x .* (1 - x);
  dy = (f / (2*T)) * (1 - 2 .*x);
endfunction

function [y, dy] = exactResult2(k, x)
  num = exp((x - 1) ./ k);
  denom = (exp(-1/k) - 1);
%  y = (num - 1) ./ denom;
  
  t1 = exp(x./k);
  t2 = exp(1/k);
  y = (t1 - t2) / (1 - t2);
  dy = (t1) .* (1/k) ./ (1 - t2);
  
  if (isnan(y))
    x
    k
    (x-1) / k
    exp((x-1)/k)
    num
    denom
    y
  endif
%  y = (num ./ denom) - 1/denom;
%  dy = (1 / (k * denom)) * num;
endfunction

function [y, dy] = exactResultAdjoint(k, x_0, x)
  
  C0 = (exp(-x_0 / k) - exp(-1/k)) / (exp(-1/k) - 1);
  C1 = -C0;
  
  for i = 1:numel(x)
  
    y(i) = C0 + C1 * exp(x(i) / k) + (heaviside(x(i) - x_0) * (1 - exp((x(i) - x_0)/k)));
    dy(i) = (1/k)*C1 * exp(x(i) / k) + (-1/k) * exp((x(i) - x_0)/k) * heaviside(x(i) - x_0);
  
  endfor
  
endfunction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k, c, b and f are the parameters of the differential equation
%   ku'' + cu' + bu + f = 0
%

function [L2Norm, H1Norm, ElemL2Error, ElemH1Error] = errorNorm(nEls,nodeCoords,connect,elDof,dFreedom,pDeg,pType,u,k,c,b,f)
  %%%%%%%%%%%%%
  % Norm of error is Sum_elem ( int((uexact - uh)^2) )
  %
  % With gaussian quadrature: Sum_elem( Sum_gpoints ( weight * (uexact - uh)^2 * h/2 ) )

  %figure(2)
  x = 0:0.001:1;
  [y, dy] = exactResultAdjoint(k(x), 0.25, x);
  
%  y
  
  subplot(2,2,3)
  hold on
  axis square
  plot(x,y)
  
  subplot(2,2,4)
  hold on
  axis square
  plot(x,dy)
  
  for iElem = 1:nEls
    for i = 1:elDof(iElem)
        c = dFreedom(iElem,i);
        uL(iElem,i) = u(c);
    end
  end
  
  ElemL2Error = zeros(nEls,1);
  ElemH1Error = zeros(nEls,1);
  TotalError = 0;
  TotalDerivError = 0;
  [xiQ,wQ] = gQuad;

  for iElem = 1 : nEls
    x1 = nodeCoords(connect(iElem, 1));
    x2 = nodeCoords(connect(iElem, 2));
    h = x2 - x1;
    jac = 2 / h;
    
    nQ = 4; % order of gaussian quadrature
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
      [uexact, duexact] = exactResultAdjoint(k(xq), 0.25, xq);
      
      % Value at this gaussian point
      Error = wq * (uexact - uh)^2 * h / 2;
      DerivError = wq * (duexact - duh)^2 * h / 2;
      TotalError = TotalError + Error;
      TotalDerivError = TotalDerivError + DerivError;
      
      ElemL2Error(iElem) = ElemL2Error(iElem) + Error;
      ElemH1Error(iElem) = ElemH1Error(iElem) + Error + DerivError;
      
    end
    
    
  end

  L2Norm = sqrt(TotalError);
  L2DerivNorm = sqrt(TotalDerivError);
  H1Norm = sqrt(TotalError + TotalDerivError);

endfunction