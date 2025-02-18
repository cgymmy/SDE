%test the P1 element for advection equation
addpath('../TimeIntegrator');
addpath('../SpacialOp');
Nx=80;%number of elements
h=2*pi/Nx;
elementnodes=[0:h:2*pi];
init=zeros(2,Nx);
for i=1:Nx
    left=initval(elementnodes(i));
    right=initval(elementnodes(i+1));
    init(1,i)=(left+right)/2;
    init(2,i)=right-left;
end
T=1;
Nt=1000;
dt=T/Nt;
res=TVDRK3(init,dt,T,@Burgersflux,elementnodes);
y=zeros(1,Nx+1);
for i=1:Nx
    y(i)=res(1,i)-0.5*res(2,i);
end
y(Nx+1)=y(1);
plot(elementnodes,y);
hold on;
plot(elementnodes,initval(elementnodes));