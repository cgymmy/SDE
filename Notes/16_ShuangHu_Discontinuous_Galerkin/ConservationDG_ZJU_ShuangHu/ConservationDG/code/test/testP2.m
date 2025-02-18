%test the P2 element for advection equation
addpath('../TimeIntegrator');
addpath('../SpacialOp');
Nx=80;%number of elements
h=2*pi/Nx;
elementnodes=[0:h:2*pi];
interpnodes=[0:h/2:2*pi];
interpval=initval(interpnodes(1:2*Nx+1));
init=zeros(3,Nx);
for i=1:Nx
    x1=interpval(2*i-1);
    x2=interpval(2*i);
    x3=interpval(2*i+1);
    init(1,i)=(x1+x3)/6+2/3*x2;
    init(2,i)=x3-x1;
    init(3,i)=(x1+x3)/3-2/3*x2;
end
T=[0,0.5,1,1.5];
Nt=800;
dt=T./Nt;
for j=1:size(T,2)
    res=TVDRK3(init,dt(j),T(j),@Burgersflux,elementnodes);
    y=zeros(1,2*Nx);
    x=zeros(1,2*Nx);
    %res=postprocess(res,elementnodes);
    for i=1:Nx
        x(2*i-1)=elementnodes(i);
        x(2*i)=elementnodes(i+1);
        %[left,right]=limiter(res,i,elementnodes);
        %y(2*i-1)=right;
        %y(2*i)=left;
        y(2*i-1)=res(1,i)-0.5*res(2,i)+res(3,i);
        y(2*i)=res(1,i)+0.5*res(2,i)+res(3,i);
    end
    plot(x,y);
    hold on;
    if j==4
        plotBurgers(1.49);
    else
        plotBurgers(T(j));
    end
    hold off;
    print("../../report/figure/Burgers("+num2str(j)+").eps",'-depsc');
end
hold on;