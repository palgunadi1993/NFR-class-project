clear 
close all
clc

load Simple_mesh.mat

num_nodes = msh.nbNod;
matrix_pos1 = msh.POS(:,1:2);
matrix_nodes1 = msh.TRIANGLES(:,1:3);
num_matrix = length(matrix_nodes1);

[A,Index] = sort(sum(matrix_pos1,2));
numindex = linspace(1,length(Index),length(Index));
matrix_pos = matrix_pos1(Index,:);

matrix_nodes = matrix_nodes1;

for i=1:length(Index)
    a = find(matrix_nodes1(:,1)==Index(i,1));
    b = find(matrix_nodes1(:,2)==Index(i,1));
    c = find(matrix_nodes1(:,3)==Index(i,1));
    matrix_nodes(a,1) = i;
    matrix_nodes(b,2) = i;
    matrix_nodes(c,3) = i;
end

sumnodes = sum(matrix_nodes,2);
[A,I] = sort(sumnodes);
matrix_nodes = matrix_nodes(I,:);

% Finding the center of the triangles
[center_matrix, matrix_id, matrix_area] = centerNodes(num_matrix, matrix_pos, matrix_nodes);

% Finding the center of connectivities and 1D-area
[center_fracture, fracture_id, frac_area, fracture_nodes_pairs] = centerLines(num_matrix, matrix_nodes, matrix_pos, matrix_id);

% Neighboring points -- matrix-matrix -> not necessary?
matrix_fracture_ID_test = 31;
[neigh_matrix_matrix_ID] = neighboringElements(matrix_nodes, matrix_fracture_ID_test);

% Matrix-fracture
[neigh_matrix_fracture_ID] = neighboringElementFracture(matrix_nodes,fracture_nodes_pairs, matrix_fracture_ID_test);

% Fracture-matrix-fracture
ID_fracture_test = num_matrix + 3;
[neigh_fracture_fracture_ID, neigh_fracture_matrix_ID] = fracture_matrix_fracture(ID_fracture_test, fracture_id, fracture_nodes_pairs, matrix_nodes);

%% to see the fracture id
figure(1)
T = delaunay(matrix_pos(:,1),matrix_pos(:,2));
% plot(center_matrix(:,1), center_matrix(:,2),'k.')
hold on
triplot(T,matrix_pos(:,1),matrix_pos(:,2))
% text(center_matrix(:,1), center_matrix(:,2),num2str(matrix_id))
% plot(matrix_pos(:,1),matrix_pos(:,2),'ok','markerfacecolor','k')
% text(matrix_pos(:,1),matrix_pos(:,2),num2str(numindex'))
plot(center_fracture(:,1), center_fracture(:,2),'r.')
text(center_fracture(:,1), center_fracture(:,2),num2str(fracture_id))

%%

% Find fracture
% x_bor = 2.0;
% y_bor = [2.0,8.0];
% ab = 1;
% a = zeros(200,1);
% for i=1:length(center_fracture)
%     if center_fracture(i,1) <= x_bor + 0.025 && center_fracture(i,1) >= x_bor - 0.025
%         if center_fracture(i,2) <= y_bor(2) && center_fracture(i,2) >= y_bor(1)
%             a(ab) = fracture_id(i);
%             ab = ab + 1;
%         end
%     end
% end
% a(a==0) = [];
% y_bor = 8.0;
% x_bor = [2.5,8.0];
% ab = 1;
% b = zeros(200,1);
% for i=1:length(center_fracture)
%     if center_fracture(i,1) <= x_bor(2) && center_fracture(i,1) >= x_bor(1)
%         if center_fracture(i,2) <= y_bor + 0.025 && center_fracture(i,2) >= y_bor - 0.025
%             b(ab) = fracture_id(i);
%             ab = ab + 1;
%         end
%     end
% end
% b(b==0) = [];
% c = [a;b];

x_bor = 5.0;
y_bor = [2.0,8.0];
ab = 1;
a = zeros(200,1);
for i=1:length(center_fracture)
    if center_fracture(i,1) <= x_bor + 0.025 && center_fracture(i,1) >= x_bor - 0.025
        if center_fracture(i,2) <= y_bor(2) && center_fracture(i,2) >= y_bor(1)
            a(ab) = fracture_id(i);
            ab = ab + 1;
        end
    end
end
a(a==0) = [];

% Initiation of fracture
% frac = [161, 182, 201, 223, 251, 278, 299, 314];
% frac = [210, 231,250,273,306,325,344,188,214,244, 272,304,341,362,380];
% frac = [1148,1100,1050,1006,974,948,931,924,900,916,952,971,993,1025,1064,1119,1164,1215,1263,...
%     1282,1315,1321,1310,1300,1156,1083,1040,982,929,876,835,796,764,741];

% frac = [411,444,490,526,559];
frac = a;% cmg

init_fracture = frac - num_matrix;
x_frac_init = matrix_pos(fracture_nodes_pairs(init_fracture,:),1);
y_frac_init = matrix_pos(fracture_nodes_pairs(init_fracture,:),2);

% change fracture_id
new_fracture_nodes_pairs = zeros(length(init_fracture),2);
new_center_fracture = zeros(length(init_fracture),2);
new_fracture_id = zeros(length(init_fracture),1);

new_fracture_id(:,1) = linspace(1+num_matrix, length(init_fracture)+num_matrix,length(init_fracture));
new_center_fracture(:,1) = center_fracture(init_fracture,1);
new_center_fracture(:,2) = center_fracture(init_fracture,2);
new_frac_area = frac_area(init_fracture,1);
new_fracture_nodes_pairs(:,1) = fracture_nodes_pairs(init_fracture,1);
new_fracture_nodes_pairs(:,2) = fracture_nodes_pairs(init_fracture,2);

% Merge_matrix_fracture id
matrix_frac_id = [matrix_id;new_fracture_id];
matrix_frac_center = [center_matrix;new_center_fracture];
matrix_frac_area = [matrix_area; new_frac_area];
con_frac = [new_fracture_nodes_pairs, zeros(length(new_fracture_nodes_pairs),1)];
matrix_frac_nodes = [matrix_nodes;con_frac];

% Fracture-matrix-fracture
ID_fracture_test = new_fracture_id(2);
[neigh_fracture_fracture_ID, neigh_fracture_matrix_ID] = fracture_matrix_fracture(ID_fracture_test, new_fracture_id, new_fracture_nodes_pairs, matrix_nodes);

% Matrix-matrix-fracture
matrix_fracture_ID_test = 69;
[neigh_matrix_matrix_ID] = neighboringElementsFrac(matrix_nodes, new_fracture_nodes_pairs, num_matrix, matrix_fracture_ID_test);


T = delaunay(matrix_pos(:,1),matrix_pos(:,2));

% for i=1:length(center_nodes(:,1))
figure(2)
plot(center_matrix(:,1), center_matrix(:,2),'k.')
hold on
triplot(T,matrix_pos(:,1),matrix_pos(:,2))
text(center_matrix(:,1), center_matrix(:,2),num2str(matrix_id))
plot(new_center_fracture(:,1), new_center_fracture(:,2),'r.')
text(new_center_fracture(:,1), new_center_fracture(:,2),num2str(new_fracture_id))

% fracture-position
plot(x_frac_init, y_frac_init,'r', 'LineWidth',2)

% % matrix-fracture
plot(center_matrix(matrix_fracture_ID_test,1), center_matrix(matrix_fracture_ID_test,2),'square')
% plot(center_fracture(neigh_matrix_fracture_ID,1), center_fracture(neigh_matrix_fracture_ID,2),'k^')

% plot(matrix_pos(matrix_nodes(matrix_fracture_ID_test,:),1), matrix_pos(matrix_nodes(matrix_fracture_ID_test,:),2),'o')
plot(matrix_frac_center(neigh_matrix_matrix_ID,1), matrix_frac_center(neigh_matrix_matrix_ID,2),'ko')

% % fracture-matrix-fracture
plot(matrix_frac_center(ID_fracture_test,1), matrix_frac_center(ID_fracture_test,2),'sk')
plot(matrix_frac_center(neigh_fracture_matrix_ID,1), matrix_frac_center(neigh_fracture_matrix_ID,2),'ko')
plot(matrix_frac_center(neigh_fracture_fracture_ID,1), matrix_frac_center(neigh_fracture_fracture_ID,2),'^k')

% hold off
% pause(1)
% end