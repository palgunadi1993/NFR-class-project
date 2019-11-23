% ******************************************** Calculate  porduction well mobilities  ****************************************************
function  [WellMobW,WellMobO]=CalculateWellMob(matrix_frac_center,num_matrix,Swij,swir,wexp,krwro,sorw,oexp,krocw,mu_o,mu_w,Wells)

N = length(matrix_frac_center);
for i=1:N
    if(Wells(i,1).id==-1)
        sw_up=Swij(i);
        [krw,kro]=CalculateRelPerm(i,i,num_matrix,sw_up,swir,sorw,wexp,oexp,krwro,krocw);
        WellMobO(i,1)=kro/mu_o;
        WellMobW(i,1)=krw/mu_w;
    end
end
