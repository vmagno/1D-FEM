clc
clear all
close all

degree = 1;
viscosity = 0.1;

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

err = [];
dofs = [];

for i=1:20
	residualError = main_residual(mMesh,degree,viscosity);
	err = [err sqrt(sum(residualError.^2))];
	dofs = [dofs mMesh.nN];
	maxErr = max(residualError);
    mMesh = refineMesh(mMesh,find(residualError>0.5*maxErr));
end

figure;
plot(log(dofs), log(err));
