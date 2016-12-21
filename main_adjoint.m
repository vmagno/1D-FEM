function resid = main_adjoint(mesh, degree, viscosity)
%*******************main code for the 1D finite element solver*************

%**************************build mesh**************************************
x_0 = 0.25;

nN = mesh.nN;
nEls = mesh.nEls;
nodeCoords = mesh.nodeCoords;
connect = mesh.connect;
nB = mesh.nB;
bEls = mesh.bEls;
bPts = mesh.bPts;

figure(2)
a = [0:1/(nN-1):1];
a = zeros(nN);
plot(nodeCoords(:,1), a, 'b-o')

%*******************degrees of freedom*************************************

pDeg=zeros(nEls,1);
pDeg(:,1)=degree;     %set polynomial degree for each element
%pDeg=[1,2,2,1]';   %could set each element individually
pType=zeros(nEls,1);
pType(:,1)=2;   %set element type: 1=Lagrangian, 2=hierarchical

[elDof,dFreedom]=dof(nEls,pDeg,connect);

%*****************build element matrices k and vectors f*******************

%get quadrature points and weights
[xiQ,wQ]=gQuad;

%build k and f for each element; sum to find K and F
%set k,c,b,f:
syms x
k=@(x)-viscosity;
c=@(x)1;
b=@(x)0;
f=@(x)-1;
%f=@(x)(heaviside(x-.3)-heaviside(x-.6))*-1;
[K,F]=element(nEls,nodeCoords,connect,xiQ,wQ,pDeg,pType,elDof,dFreedom,k,c,b,f);

%*****************add boundary conditions and solve************************

%set boundary condition type (set values within boundaryC.m)
bType=zeros(2,1);   %boundary condition types; use 1 for Dirichlet,
bType(1)=1;         %2 for Neumann, 3 for Robin
bType(2)=1;

[K,F]=boundaryC(nB,bEls,bPts,bType,dFreedom,K,F);

%F = addDirac(F, nodeCoords, connect, x_0, pDeg, pType, dFreedom, elDof);

u=K\F;

%**************post-processing*********************************************
%plot the solution u and its derivative
figure(1)
postProc(nEls,nodeCoords,connect,elDof,dFreedom,pDeg,pType,u);
%[L2Norm, H1Norm, ElemL2Error, ElemH1Error] = errorNorm(nEls,nodeCoords,connect,elDof,dFreedom,pDeg,pType,u,k,c,b,f,x_0);
adj = @(x)adjoint_new(viscosity,x_0,x);
resid = residual_compute(nEls,nodeCoords,connect,elDof,dFreedom,pDeg,pType,u,adj,k,c,b,f);
end
