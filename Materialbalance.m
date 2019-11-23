function     [wip,cumwi,cumwp,err]=Materialbalance (matrix_frac_center,Swij,Wells,Prodw,PorVol,dt,cumwi,cumwp,owip);

N = length(matrix_frac_center);

wip=0;

for i=1:N
    wip=wip+PorVol(i,1)*Swij(i,1);
    cumwi=cumwi+Wells(i,1).rate*dt;
    cumwp=cumwp+Prodw(i,1)*dt;
end
err=abs((wip-owip)-(cumwi+cumwp));