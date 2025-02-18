%test the P2 element for advection equation
addpath('../TimeIntegrator');
addpath('../SpacialOp');
Nx=[20,40,80,160,320];%number of elements
T=0.2;
h=2*pi./Nx;
dt=0.1*h;
n=size(Nx,2);
err=zeros(1,n);
rate=zeros(1,n-1);
for i=1:n
    elementnodes=0:h(i):2*pi;
    init=getinitval(elementnodes,2);
    res=TVDRK3(init,dt(i),T,@Burgersflux,elementnodes);
    err(i)=Linferr(res,@exactBurgers,elementnodes,T);
end
for i=1:n-1
    rate(i)=log(err(i)/err(i+1))/log(2);
end
err
rate