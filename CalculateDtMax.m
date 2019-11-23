function [dtnew,icut]=DtControl(Nx,Ny,dt,Swij_n,Swij,dtmin,dtmax,dswmin,dswmax);
icut=0;
dtnew=dt;
dsw=0;
for j=1:Ny
    for i=1:Nx
        dsw=max(abs(Swij_n(i,j)-Swij(i,j)),dsw);
    end
end

if(dsw>dswmax) 
    dtnew=dt*0.5;
    dtnew=max(dtmin,dtnew);
    if (dtnew>dtmin)icut=1;end
elseif(dsw<dsmin)
    dtnew=dt*1.02;
    dtnew=min(dtmax,dtnew);
end