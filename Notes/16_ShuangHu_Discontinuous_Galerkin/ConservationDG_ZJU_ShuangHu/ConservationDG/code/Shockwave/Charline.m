addpath('../test');
Lines=10;
h=2*pi/Lines;
xzero=0:h:2*pi;
t=0:0.1:1.2;
for i=1:size(xzero,2)
    start=xzero(i);
    k=initval(start);
    x=start+k*t;
    plot(x,t);
    hold on;
end
print('../../report/figure/Charline.eps','-depsc');
