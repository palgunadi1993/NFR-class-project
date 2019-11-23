function [center_nodes, num_id, area] = centerNodes(num_triangles, pos_nodes, tri_nodes)
% Finding the center of triangles

center_nodes = zeros(num_triangles,2);
num_id = zeros(num_triangles,1);
area = zeros(num_triangles,1);
for i=1:num_triangles
    center_nodes(i,1) = mean(pos_nodes(tri_nodes(i,:),1));
    center_nodes(i,2) = mean(pos_nodes(tri_nodes(i,:),2));
    tempx = pos_nodes(tri_nodes(i,:),1);
    tempy = pos_nodes(tri_nodes(i,:),2);
    a = abs(pdist([tempx(1), tempy(1);tempx(2), tempy(2)],'euclidean'));
    b = abs(pdist([tempx(1), tempy(1);tempx(3), tempy(3)],'euclidean'));
    c = abs(pdist([tempx(2), tempy(2);tempx(3), tempy(3)],'euclidean'));
    s = (a+b+c)/2;
    area(i) = sqrt(s.*(s-a).*(s-b).*(s-c));
    num_id(i,1) = i;
end