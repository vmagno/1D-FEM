%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% k, c, b and f are the parameters of the differential equation
%   ku'' + cu' + bu + f = 0
%
function res = residual_compute(nEls,nodeCoords,connect,elDof,dFreedom,pDeg,pType,u,adjoint,k,c,b,f)
  % Returns residual per element

  for iElem = 1:nEls
    for i = 1:elDof(iElem)
        glob = dFreedom(iElem,i);
        uL(iElem,i) = u(glob);
    end
  end
  
  res = [];
  nQ = 3; % order of gaussian quadrature
  [xiQ,wQ] = gQuad; % quadrature points and weights

  for iElem = 1 : nEls
	res_local = 0;

    x1 = nodeCoords(connect(iElem, 1));
    x2 = nodeCoords(connect(iElem, 2));
	h = x2 - x1;
    jac = 2/h;

	% Compute L2 norm of element residual using Gaussian quadrature
    for iG = 1:nQ;
      xi = xiQ(iG,nQ);
      wq = wQ(iG,nQ);
      
      % Get value of computed solution at point xi
      [N,dN] = shape(xi,iElem,pDeg,pType);
	  dN = dN * jac;
	  uh = 0;
      duh = 0;
      for iDof = 1:elDof(iElem)
        duh = duh + dN(iDof) * uL(iElem, iDof);
      end
	  xi_glob = xi*h/2 +(x1+x2)/2;
	  [p,dP] = adjoint(xi_glob);
      
	  % Residual value at xi
	  rk = p -k(xi)*duh*dP - c(xi)*duh*p;
      % Gaussian cubature point for L2 norm of rk
	  res_local += wq * rk;
    end

	res = [res, res_local];
  end
end
