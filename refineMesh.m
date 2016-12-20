function outMesh = refineMesh(mesh, elementIds)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refine the input mesh by dividing the chosen elements into 2 elements each
%
% Takes the initial number of nodes, their coordinates, the connectivity and
% an array of elements to divide
%
% Returns the final number of nodes
%

outMesh = mesh;

for iElem = unique(elementIds)
  
  % Midpoint
  elemMiddle = (mesh.nodeCoords(mesh.connect(iElem,1)) + mesh.nodeCoords(mesh.connect(iElem,2))) / 2;
  
  % Add the coordinates
  outMesh.nodeCoords(++outMesh.nN) = elemMiddle;

  % Adjust connectivity
  outMesh.connect(iElem,2) = outMesh.nN;
  outMesh.connect = [outMesh.connect; [outMesh.nN, mesh.connect(iElem,2)]];
  
  outMesh.nEls++;

  bIdx = find(mesh.bEls==iElem);
  if bIdx & mesh.bPts(bIdx)==2
	  outMesh.bEls(bIdx)=outMesh.nEls;
  end
end

end
