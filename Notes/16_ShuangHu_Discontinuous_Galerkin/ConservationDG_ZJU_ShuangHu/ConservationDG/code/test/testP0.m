%test the P0 element for advection equation
addpath('../TimeIntegrator');
addpath('../SpacialOp');
Nx=800;%number of elements
h=2*pi/Nx;
elementnodes=[0:h:2*pi];
init=initval(elementnodes(1:Nx));
T=2*pi;
Nt=1000;
dt=T/Nt;
res=TVDRK3(init,dt,T,@advectionflux,elementnodes);
res(1,Nx+1)=res(1);
plot(elementnodes,res);
hold on;
plot(elementnodes,initval(elementnodes));