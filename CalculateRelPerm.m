
% -----------------------------------------------------------
function [krw,kro]=CalculateRelPerm(center_id,neigh_ID,num_matrix,sw,swir,sorw,wexp,oexp,krwro,krocw) 

if(sw<=swir)
    krw=0;
    kro=krocw;
    return
end
if(sw>=(1-sorw))
    krw=krwro;
    kro=0;
    return; 
end

if center_id > num_matrix && neigh_ID > num_matrix
    % Fracture
    krw = 2.22222*sw - 0.55555;
    kro = -2.22222*sw + 1.55554;
else
%     matrix
    krw = krwro * ((sw - swir) / (1 - swir - sorw))^wexp;
    kro = krocw * ((1- sw - sorw) / (1 - swir - sorw))^oexp;
end
