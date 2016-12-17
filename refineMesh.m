function [outNumNodes, outNumElem, outNodeCoords, outConnect] = refineMesh(nNodes, nodeCoords, connect, lastElem, elementIds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine the input mesh by dividing the chosen elements into 2 elements each
%
% Takes the initial number of nodes, their coordinates, the connectivity and
% an array of elements to divide
%
% Returns the final number of nodes
%

outNumNodes = nNodes;
outNodeCoords = nodeCoords;
outConnect = connect;

for iToDivide = 1:numel(elementIds)
  iElem = elementIds(iToDivide);
  
  % Midpoint
  elemMiddle = (nodeCoords(connect(iElem,1)) + nodeCoords(connect(iElem,2))) / 2;
  
  % Add the coordinates
  outNumNodes = outNumNodes + 1;
  outNodeCoords(outNumNodes) = elemMiddle;
  
  % Adjust connectivity
  if (iElem != lastElem)
    outConnect(outNumNodes-1,2) = connect(iElem,2);
    outConnect(iElem,2) = outNumNodes;
    outConnect(outNumNodes-1,1) = outNumNodes;
  else
    outConnect(outNumNodes-1,1) = connect(iElem,1);
    outConnect(iElem,1) = outNumNodes;
    outConnect(outNumNodes-1,2) = outNumNodes;
  endif
  
  outNodeCoords
  outConnect
  
endfor

outNumElem = outNumNodes - 1;

endfunction