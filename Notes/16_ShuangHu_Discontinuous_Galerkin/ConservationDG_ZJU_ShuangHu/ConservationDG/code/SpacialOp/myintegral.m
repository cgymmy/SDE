function res=myintegral(f,left,right)
    res=0;
    coe=(right-left)/2;
    weight=[0.17132449,0.36076157,0.46791393,0.46791393,0.36076157,0.17132449];
    point=[-0.93246951,-0.66120929,-0.23861919,0.23861919,0.66120929,0.93246951];
    point=(left+right)/2+(right-left)/2*point;
    for i=1:6
        res=res+weight(i)*f(point(i));
    end
    res=res*coe;
end