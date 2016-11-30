function postProc(nEls,nodeCoords,connect,elDof,dFreedom,pDeg,pType,u)
%plot the solution u and its derivative

for ne = 1:nEls
    for i = 1:elDof(ne)
        c = dFreedom(ne,i);
        uL(ne,i) = u(c);
    end
end
for ne = 1:nEls
    nx = 10;
    u = zeros(nx,1);
    du = zeros(nx,1);
    x1 = nodeCoords(connect(ne,1));
    x2 = nodeCoords(connect(ne,2));
    h = x2 - x1;
    jac = 2/h;
    x = (x1:(x2-x1)/(nx-1):x2);
    xi = (-1:2/(nx-1):1);
    for j = 1:nx
        [N,dN] = shape(xi(j),ne,pDeg,pType);
        dN = dN*jac;
        for i = 1:elDof(ne)
            u(j) = u(j) + N(i)*uL(ne,i);
            du(j) = du(j) + dN(i)*uL(ne,i);
        end
    end
    
    subplot(2,2,1)
    hold on
    axis square
    plot(x,u)
    subplot(2,2,2)
    hold on
    axis square
    plot(x,du)
    
    subplot(2,2,3)
    hold on
    axis square
    plot(x,exactResult(-10, 10, x))
    
end

%%%%%%%%%%%%%
% Norm of error is Sum_elem ( int((uexact - uh)^2) )
%
% With gaussian quadrature: Sum_elem( Sum_gpoints ( weight * (uexact - uh)^2 * h/2 ) )

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
    [N,dN] = shape(xi,ne,pDeg,pType);
    uh = 0;
    duh = 0;
    for iDof = 1:elDof(iElem)
      uh = uh + N(iDof) * uL(iElem, iDof);
      duh = duh + dN(iDof) * uL(iElem, iDof) * 2 / h;
    end
    
    % Value of exact solution at point xi (xq)
    [uexact, duexact] = exactResult(-10, 10, xq);
    
    % Value at this gaussian point
    TotalError = TotalError + wq * (uexact - uh)^2 * h / 2;
    TotalDerivError = TotalDerivError + wq * (duexact - duh)^2 * h / 2;
    
  end
  
  
end

L2Norm = sqrt(TotalError)
L2DerivNorm = sqrt(TotalDerivError)