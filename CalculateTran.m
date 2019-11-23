% ******************************************** Calculate  transmissibilities ****************************************************
function [Trans]=CalculateTran(matrix_frac_center,matrix_nodes,matrix_frac_nodes, new_fracture_nodes_pairs, num_matrix,new_fracture_id,matrix_pos,aperture,Permx, max_intersection)
% Transmissibilities 

N = length(matrix_frac_center);
Trans=zeros(N,4+max_intersection);

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
        frac_neigh = neigh_ID(neigh_ID>num_matrix);
        nodes_pt = matrix_frac_nodes(i,:);
        temp1 = matrix_frac_nodes(frac_neigh,:);
        temp1 = temp1(:,1:2);
        s1 = sum(sum(ismember(temp1,nodes_pt(1,1)),2));
        if s1 > 1
            temp2 = nodes_pt(1,1);
        else
            temp2 = nodes_pt(1,2);
        end
        new_frac_neigh = frac_neigh(logical(sum(ismember(temp1,temp2),2)),:);
        Dwf = 0;
        for l=1:length(new_frac_neigh)
            pt_ctr = matrix_frac_center(i,:);
            nodes = matrix_frac_nodes(new_frac_neigh(l),:);
            con_pt = nodes_pt(ismember(nodes_pt,nodes));
            con_pt(con_pt==0) = [];
            
            D1 = abs(pdist([pt_ctr;matrix_pos(con_pt(1),:)],'euclidean'));
            D2 = abs(pdist([matrix_pos(con_pt(1),:);matrix_frac_center(new_frac_neigh(l),:)],'euclidean'));
            
            if l==1
                t1=Permx(i,1)/(D1);
            else
                t1=0;
            end
            t2=Permx(new_frac_neigh(l),1)/(D2);
            Dwf = Dwf + t1 + t2;
        end
    end
    for j=1:num_neigh
        % normal distance initiation
        nodes_pt = matrix_frac_nodes(i,:);
        pt_ctr = matrix_frac_center(i,:);
        nodes = matrix_frac_nodes(neigh_ID(j),:);
        con_pt = nodes_pt(ismember(nodes_pt,nodes));
        con_pt(con_pt==0) = [];
        
        % if the connection is not fracture-fracture
        if length(con_pt) > 1
            v1 = matrix_pos(con_pt(1),:);
            v2 = matrix_pos(con_pt(2),:);
            nei_ctr = matrix_frac_center(neigh_ID(j),:);

            % center node
            D1 = point_to_line([pt_ctr,0],[v1,0],[v2,0]);
            % neighboring nodes
            D2 = point_to_line([nei_ctr,0],[v1,0],[v2,0]);

            % if neighboring point is fracture and center node is matrix
            if (neigh_ID(j) > num_matrix)
                D2 = aperture*0.5;
            % if the center node is fracture
            elseif (i > num_matrix)
                D1 = aperture*0.5;
            end
            A = abs(pdist([v1;v2],'euclidean'));
        
            t1=Permx(i,1)/(D1);
            t2=Permx(neigh_ID(j),1)/(D2);
            U = t1*t2;
            Dw = t1 + t2;
            
            
        else
            % fracture-fracture connection
            D1 = abs(pdist([pt_ctr;matrix_pos(con_pt(1),:)],'euclidean'));
            D2 = abs(pdist([matrix_pos(con_pt(1),:);matrix_frac_center(neigh_ID(j),:)],'euclidean'));
            A = aperture;

            t1=Permx(i,1)/(D1);
            t2=Permx(neigh_ID(j),1)/(D2);
            U = t1*t2;
            
            % fracture at intersection
            if num_neigh > 4 && con_pt(1) == temp2
                Dw = Dwf;
            else
                Dw = t1 + t2;
            end
        end
        
        Trans(i,j)=A*U/(Dw);
    end
    if num_neigh > 4
        update_num = num_neigh + 1;
        for z=1:length(new_frac_neigh)
            c = z + 1;
            for d=c:length(new_frac_neigh)
                nodes_pt = matrix_frac_nodes(new_frac_neigh(z),:);
                pt_ctr = matrix_frac_center(new_frac_neigh(z),:);
                nodes = matrix_frac_nodes(new_frac_neigh(d),:);
                con_pt = nodes_pt(ismember(nodes_pt,nodes));
                con_pt(con_pt==0) = [];
                
                D1 = abs(pdist([pt_ctr;matrix_pos(con_pt(1),:)],'euclidean'));
                D2 = abs(pdist([matrix_pos(con_pt(1),:);matrix_frac_center(new_frac_neigh(d),:)],'euclidean'));
                
                A = aperture;

                t1=Permx(new_frac_neigh(z),1)/(D1);
                t2=Permx(new_frac_neigh(d),1)/(D2);
                U = t1*t2;
                
                Trans(i,update_num) = A * U/Dwf;
                update_num = update_num + 1;
            end
        end
    end
end
