%get the numerical flux
function numflux=NumericalFlux(flux,u,left,right,elementid,x)
    ux=assembleu(u,x,left,right,elementid);
    numflux=flux(ux);
    return;
end