% ******************************************** Calculate  transmissibilities ****************************************************
function  [Fluxw,Prodw,ProdRates]=CalculateWflux(matrix_frac_center,matrix_frac_nodes,matrix_nodes,new_fracture_nodes_pairs,num_matrix,new_fracture_id,Pij,Trans,Swij,swir,wexp,krwro,sorw,oexp,krocw,mu_o,mu_w,Wells,PI,Fluxw,Prodw,ProdRates);

%Fluxw_x(1:Nx+1,1:Ny)=0;
%Fluxw_y(1:Nx,1:Ny+1)=0;
%Prodw(1:Nx,1:Ny)=0;
N = length(matrix_frac_center);
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
        new_frac_neigh = frac_neigh(logical(sum(ismember(temp1,temp2),2)),:);
    end
    
    for j=1:num_neigh
        p1=Pij(i);
        p2=Pij(neigh_ID(j));
        if(p2-p1>0)
            sw_up=Swij(neigh_ID(j));
        else
            sw_up=Swij(i);
        end
        [krw,kro]=CalculateRelPerm(i,neigh_ID(j),num_matrix,sw_up,swir,sorw,wexp,oexp,krwro,krocw);
        lambda_o=kro/mu_o *0;
        lambda_w=krw/mu_w;
        Fluxw(i,j)=Trans(i,j)*(lambda_o+lambda_w)*(p2-p1);
    end
    
    if num_neigh > 4
        update_num = num_neigh + 1;
        for z=1:length(new_frac_neigh)
            c = z + 1;
            for d=c:length(new_frac_neigh)
                p1=Pij(new_frac_neigh(z));
                p2=Pij(new_frac_neigh(d));
                if(p2-p1>0)
                    sw_up=Swij(new_frac_neigh(d));
                else
                    sw_up=Swij(new_frac_neigh(z));
                end
                [krw,kro]=CalculateRelPerm(new_frac_neigh(z),new_frac_neigh(d),num_matrix,sw_up,swir,sorw,wexp,oexp,krwro,krocw);
                lambda_o=kro/mu_o *0;
                lambda_w=krw/mu_w;
                Fluxw(i,update_num)=Trans(i,update_num)*(lambda_o+lambda_w)*(p2-p1);
                update_num = update_num + 1;
            end
        end
    end
    
end



% production well
ProdRates(1:2)=0;
for i=1:N
    if (Wells(i,1).id==-1)
        sw_up=Swij(i,1);

        p1=Pij(i,1);
        p2=Wells(i,1).bhp;
        [krw,kro]=CalculateRelPerm(i,i,num_matrix,sw_up,swir,sorw,wexp,oexp,krwro,krocw);
        lambda_o=kro/mu_o;
        lambda_w=krw/mu_w;
        Prodw(i,1)=PI*(lambda_w)*(p2-p1);
        ProdRates(1)=ProdRates(1)-PI*(lambda_w)*(p2-p1);
        ProdRates(2)=ProdRates(2)-PI*(lambda_o)*(p2-p1);
    end
end



