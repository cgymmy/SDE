function numflux=LaxFriedrichesFlux(flux,u,eid,elementnodes)
    alpha=1.5;
    uminus=assembleu(u,elementnodes(eid+1),elementnodes(eid),elementnodes(eid+1),eid);
    n=size(elementnodes,2);
    if eid+1==n
        newid=1;
    else
        newid=eid+1;
    end
    uplus=assembleu(u,elementnodes(newid),elementnodes(newid),elementnodes(newid+1),newid);
    numflux=0.5*(flux(uminus)+flux(uplus)-alpha*(uplus-uminus));
    return;
end