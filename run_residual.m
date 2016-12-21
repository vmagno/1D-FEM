clc
clear all
close all

viscosity = 0.0001
for degree = 1:2
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

	err{degree} = [];
	dofs{degree} = [];

	for i=1:20
		residualError = main_residual(mMesh,degree,viscosity);
		err{degree} = [err{degree} sqrt(sum(residualError.^2))];
		dofs{degree} = [dofs{degree} mMesh.nN];
		maxErr = max(residualError);
		mMesh = refineMesh(mMesh,find(residualError>0.5*maxErr));
	end
end

figure;
plot(log(dofs{1}), log(err{1}), 'b-x;éléments linéaires;', log(dofs{2}),log(err{2}),'k--o;éléments quadratiques;');
xlabel('log(DoF)');
ylabel('log(||e||_1)');
title(sprintf('epsilon = %d', viscosity));
