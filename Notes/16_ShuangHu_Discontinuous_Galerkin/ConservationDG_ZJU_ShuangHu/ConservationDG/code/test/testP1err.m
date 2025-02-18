%test the P1 element for advection equation
addpath('../TimeIntegrator');
addpath('../SpacialOp');
Nx=[8,16,32,64,128];%number of elements
T=0.1;
Nt=[200,400,800,1600,3200];
dt=T./Nt;
n=size(Nx,2);
err=zeros(1,n);
rate=zeros(1,n-1);
for i=1:n
    h=2*pi/Nx(i);
    elementnodes=0:h:2*pi;
    init=getinitval(elementnodes,1);
    res=TVDRK3(init,dt(i),T,@Burgersflux,elementnodes);
    err(i)=L1err(res,@exactBurgers,elementnodes,T);
end
for i=1:n-1
    rate(i)=log(err(i)/err(i+1))/log(2);
end
err
rate