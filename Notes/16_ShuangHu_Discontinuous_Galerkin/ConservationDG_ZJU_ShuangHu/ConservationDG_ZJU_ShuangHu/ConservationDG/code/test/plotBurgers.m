function void=plotBurgers(T)
Nx=2000;
dx=2*pi/Nx;
x=0:dx:2*pi;
for i=1:Nx+1
    res(i)=exactBurgers(x(i),T);
end
plot(x,res);
end