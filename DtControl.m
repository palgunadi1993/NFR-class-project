function [dtnew,icut]=DtControl(matrix_frac_center,dt,Swij_n,Swij,dtmin,dtmax,dswmin,dswmax);
N = length(matrix_frac_center);
icut=0;
dtnew=dt;
dsw=0;
for i=1:N
    if(dsw<abs(Swij_n(i,1)-Swij(i,1))) 
     ijloc=[i 1];
     dsw=abs(Swij_n(i,1)-Swij(i,1));
    end
end
%ijloc
%dsw


if(dsw>dswmax) 
    dtnew=dt*0.5;
    dtnew=max(dtmin,dtnew);
    if (dtnew>dtmin)icut=1;end
elseif(dsw<dswmin)
    dtnew=dt*1.05;
    dtnew=min(dtmax,dtnew);
end