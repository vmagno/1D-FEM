clc
clear all
close all

viscosity = 0.01
for degree = 2
	nEls = 2;
	[nN,nodeCoords,connect,nB,bEls,bPts]=mesh(0,1,nEls);
	mMesh = struct(
	'nN', nN,
	'nEls', nEls,
	'nodeCoords',nodeCoords,
	'connect',connect,
	'nB',nB,
	'bEls',bEls,
	'bPts',bPts);

	resid_total = [];

	for i=1:10
		resid = main_adjoint(mMesh,degree,viscosity);
		maxErr = max(abs(resid));
		resid_total = [resid_total sum(resid)];
		mMesh = refineMesh(mMesh,find(abs(resid)>0.5*maxErr));
	end
end

figure;
plot(1:10, resid_total);
