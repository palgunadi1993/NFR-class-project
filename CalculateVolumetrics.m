function [ooip,owip,totpv,PorVol]=CalculateVolumetrics(matrix_frac_aper_area,Swij,por,vol_unit);

num_id = length(matrix_frac_aper_area);
PorVol=zeros(num_id,1);
for i=1:num_id
    PorVol(i,1)= matrix_frac_aper_area(i,1)*por(i,1)/vol_unit;
end


ooip=0;
owip=0;
totpv=0;
for i=1:num_id
    ooip=ooip+PorVol(i,1)*(1-Swij(i,1));
    owip=owip+PorVol(i,1)*Swij(i,1);
    totpv=totpv+PorVol(i,1);
end