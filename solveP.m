function [P,Pij]=solveP(matrix_frac_center,Trans,num_matrix,matrix_nodes,new_fracture_nodes_pairs,new_fracture_id,Wells,PI,WellMobW,WellMobO,Jac,Rhs)

N = length(matrix_frac_center);
%Jac=sparse(Nt,Nt);
%Jac=0;
%Rhs=zeros(Nt,1);


ic=0;
for i=1:N
    ic=ic+1;
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
    Jac(ic,ic)=0;
    Rhs(ic)=0;
    % if producer (i.e. pressure BC)
    if Wells(i,1).id==-1
        Jac(ic,ic)=Jac(ic,ic)+PI*(WellMobW(i,1)+WellMobO(i,1));
        Rhs(ic,1)=Rhs(ic,1)+PI*(WellMobW(i,1)+WellMobO(i,1))*Wells(i,1).bhp;
        % else if there is an injector, or no well
    end
    Rhs(ic,1)=Rhs(ic,1)+Wells(ic,1).rate;
    for j=1:num_neigh
        %
        Jac(ic,neigh_ID(j)) = -Trans(i,j); %- Pc;
        Jac(ic,ic) = Jac(ic,ic) - Jac(ic,neigh_ID(j));
    end
end
%Jac
% condest(Jac)
%Rhs
%solve the system
% size(Jac);
% size(Rhs);
% spy(Jac);
% 
%  solve the linear system

P=Jac\Rhs;

Pij = P;