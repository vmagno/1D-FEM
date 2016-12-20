clc
clear all
close all

degree = 1;
viscosity = 1;

xMin=0;
xMax=1;
nEls=2;
[nN,nodeCoords,connect,nB,bEls,bPts]=mesh(xMin,xMax,nEls);
mMesh = struct(
	'nN', nN,
	'nEls', nEls,
	'nodeCoords',nodeCoords,
	'connect',connect,
	'nB',nB,
	'bEls',bEls,
	'bPts',bPts);

for i=1:10
	residualError = main_residual(mMesh,degree,viscosity);
	maxErr = max(residualError);
    mMesh = refineMesh(mMesh,find(residualError>0.5*maxErr));
end

