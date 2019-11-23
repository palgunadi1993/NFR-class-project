function [Swij_n]=CalculateSw(matrix_frac_center,Fluxw,Prodw,Swij,Swij_n,dt,Wells,PorVol)

N = length(matrix_frac_center);
for i=1:N
    sumqw=sum(Fluxw(i,:))+Prodw(i,1);
    Swij_n(i,1)= Swij(i,1)+dt*sumqw/PorVol(i,1) + dt*Wells(i,1).rate/PorVol(i,1);
end


% function [Swij_n]=CalculateSw(matrix_frac_center,matrix_nodes,new_fracture_nodes_pairs,num_matrix,Fluxw,new_fracture_id,Prodw,Swij,Swij_n,dt,Wells,PorVol)
% 
% N = length(matrix_frac_center);
% for i=1:N
%     % Finding neighboring points
%     if i <= num_matrix
%         % center node (i) is a matrix
%         [neigh_ID] = neighboringElementsFrac(matrix_nodes, new_fracture_nodes_pairs, num_matrix, i);
%     else
%         % center node (i) is a fracture
%         [neigh_fracture_fracture_ID, neigh_fracture_matrix_ID] = fracture_matrix_fracture(i, new_fracture_id, new_fracture_nodes_pairs, matrix_nodes);
%         neigh_ID = [neigh_fracture_fracture_ID;neigh_fracture_matrix_ID];
%     end
%     num_neigh = length(neigh_ID);
%     for j=1:num_neigh
%         sumqw=sum(Fluxw(i,:));
%     end
% %     +Prodw(i)
%     Swij_n(i,1)= Swij(i,1)+dt*sumqw/PorVol(i,1) + dt*Wells(i,1).rate/PorVol(i,1);
% end
% 
% end