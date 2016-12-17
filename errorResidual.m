%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k, c, b and f are the parameters of the differential equation
%   ku'' + cu' + bu + f = 0
%
function err = errorResidual(nEls,nodeCoords,connect,elDof,dFreedom,pDeg,pType,u,k,c,b,f)
  % Returns residual error upper-bound per element

  for iElem = 1:nEls
    for i = 1:elDof(iElem)
        glob = dFreedom(iElem,i);
        uL(iElem,i) = u(glob);
    end
  end
  
  err = [];
  nQ = 3; % order of gaussian quadrature
  [xiQ,wQ] = gQuad; % quadrature points and weights

  for iElem = 1 : nEls
	err_local = 0;

    x1 = nodeCoords(connect(iElem, 1));
    x2 = nodeCoords(connect(iElem, 2));
	h = x2 - x1;
    jac = 2/h;

	% Compute L2 norm of element residual using Gaussian quadrature
    for iG = 1:nQ;
      xi = xiQ(iG,nQ);
      wq = wQ(iG,nQ);
      
      % Get value of computed solution at point xi
      [N,dN,d2N] = shape(xi,iElem,pDeg,pType);
	  dN = dN * jac;
	  d2N = d2N * jac * jac;
	  uh = 0;
      duh = 0;
      d2uh = 0;
      for iDof = 1:elDof(iElem)
        uh = uh + N(iDof) * uL(iElem, iDof);
        duh = duh + dN(iDof) * uL(iElem, iDof);
        d2uh = d2uh + d2N(iDof) * uL(iElem, iDof);
      end
      
	  % Residual value at xi
	  rk = k(xi)*d2uh + c(xi)*duh + b(xi)*uh + f(xi);
      % Gaussian cubature point for L2 norm of rk
	  err_local = err_local + wq * rk^2;
    end

	err = [err, h^2 * err_local*(-1/k(xi))];
  end
end
