function [neighElement_id] = neighboringElements(tri_nodes, centerElement_id)
% Finding element matrix

%% Sorting similar ids
a = tri_nodes(centerElement_id,:);
ab1 = find(tri_nodes(:,1)==a(1));
ab2 = find(tri_nodes(:,2)==a(1));
ab3 = find(tri_nodes(:,3)==a(1));
ab = sort([ab1;ab2;ab3]);

ac1 = find(tri_nodes(:,1)==a(2));
ac2 = find(tri_nodes(:,2)==a(2));
ac3 = find(tri_nodes(:,3)==a(2));
ac = sort([ac1;ac2;ac3]);

ad1 = find(tri_nodes(:,1)==a(3));
ad2 = find(tri_nodes(:,2)==a(3));
ad3 = find(tri_nodes(:,3)==a(3));
ad = sort([ad1;ad2;ad3]);

%% Finding member elements
shared_elem1 = ismember(ab,ac);
shared_elem2 = ismember(ab,ad);
shared_elem3 = ismember(ac,ad);

%% Mapping neighboring elements
id_ab = ab(shared_elem1);
id_ab(id_ab==centerElement_id) = [];
id_ac = ab(shared_elem2);
id_ac(id_ac==centerElement_id) = [];
id_ad = ac(shared_elem3);
id_ad(id_ad==centerElement_id) = [];

neighElement_id = [id_ab;id_ac;id_ad];