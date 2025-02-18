function e=Linferr(res,truesol,elementnodes,T)
    n=size(elementnodes,2)-1;
    e=0;
    for i=1:n
        func=@(x)abs(assembleu(res,x,elementnodes(i),elementnodes(i+1),i)-truesol(x,T));
        e=max(e,func(elementnodes(i)));
    end
    return;
end