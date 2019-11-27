% ******************************************** Calculate  transmissibilities ****************************************************
function  [Mob]=CalculateMob(matrix_frac_center,matrix_frac_nodes,max_intersection,Pij,Trans,num_matrix,matrix_nodes,new_fracture_nodes_pairs,new_fracture_id,Mob,Swij,swir,wexp,krwro,sorw,oexp,krocw,mu_o,mu_w)

N = length(matrix_frac_center);
Mob(1:N,1:4+max_intersection)=0;
for i=1:N
    % Finding neighboring points
    if i <= num_matrix
        % center node (i) is a matrix
        [neigh_ID] = neighboringElementsFrac(matrix_nodes, new_fracture_nodes_pairs, num_matrix, i);
    else
        % center node (i) is a fracture
        [neigh_fracture_fracture_ID, neigh_fracture_matrix_ID] = fracture_matrix_fracture(i, new_fracture_id, new_fracture_nodes_pairs, matrix_nodes);
        neigh_ID = [neigh_fracture_fracture_ID;neigh_fracture_matrix_ID];
    end
    num_neigh = length(neigh_ID);
    
    % fracture at intersection
    if num_neigh > 4
        nodes_pt = matrix_frac_nodes(i,:);
        frac_neigh = neigh_ID(neigh_ID>num_matrix);
        temp1 = matrix_frac_nodes(frac_neigh,:);
        temp1 = temp1(:,1:2);
        s1 = sum(sum(ismember(temp1,nodes_pt(1,1)),2));
        if s1 > 1
            temp2 = nodes_pt(1,1);
        else
            temp2 = nodes_pt(1,2);
        end
        new_frac_neigh1 = frac_neigh(logical(sum(ismember(temp1,temp2),2)),:);
        new_frac_neigh = [i;new_frac_neigh1];
        neigh_ID(ismember(neigh_ID,new_frac_neigh)) = [];
    end
    
    for j=1:length(neigh_ID)
        p1=Pij(i); % pressure in the middle of the grid
        p2=Pij(neigh_ID(j));
        if(p2-p1>0)
            sw_up=Swij(neigh_ID(j));
        else
            sw_up=Swij(i);
        end
        [krw,kro]=CalculateRelPerm(i,neigh_ID(j),num_matrix,sw_up,swir,sorw,wexp,oexp,krwro,krocw);
        lambda_o=kro/mu_o; lambda_w=krw/mu_w;
        Mob(i,j)=Trans(i,j)*(lambda_o+lambda_w);
    end
    
    if num_neigh > 4
        update_num = length(neigh_ID) + 1;
        z = 1;
        c = z + 1;
        for d=c:length(new_frac_neigh)
            p1=Pij(new_frac_neigh(z)); % pressure in the middle of the grid
            p2=Pij(new_frac_neigh(d));
            if(p2-p1>0)
                sw_up=Swij(new_frac_neigh(d));
            else
                sw_up=Swij(new_frac_neigh(z));
            end
            [krw,kro]=CalculateRelPerm(new_frac_neigh(z),new_frac_neigh(d),num_matrix,sw_up,swir,sorw,wexp,oexp,krwro,krocw);
            lambda_o=kro/mu_o; lambda_w=krw/mu_w;
            Mob(i,update_num)=Trans(i,update_num)*(lambda_o+lambda_w);
            update_num = update_num + 1;
        end
    end
    
end

