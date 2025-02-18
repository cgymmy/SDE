function f=assembleu(u,x,left,right,elementid)
    m=size(u,1);
    f=0;
    for i=1:m
        f=f+(u(i,elementid))*basis1D(x,i-1,0,left,right);
    end
    return;
end