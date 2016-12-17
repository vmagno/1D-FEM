function outF = addDirac(F, nodeCoords, connect, x_0, pDeg, pType, dFreedom, elDof)
  
  outF = F;
  
  for iElem = 1:numel(nodeCoords)-1
    
    if (nodeCoords(connect(iElem,1)) < x_0 && nodeCoords(connect(iElem, 2)) >= x_0)
      
      [N,dN] = shape(x_0,iElem,pDeg,pType);
      
      for iDof = 1:elDof(iElem)
        iGDof = dFreedom(iElem, iDof);
        outF(iGDof) = outF(iGDof) + N(iDof);
      endfor
      
    endif
    
  endfor
  
endfunction