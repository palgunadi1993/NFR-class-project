function [neighAll_id] = neighboringElementsFrac(tri_nodes, pairs, num_matrix, centerElement_id)
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
id_ab1 = ab(shared_elem1);
id_ab1(id_ab1==centerElement_id) = [];
id_ac1 = ab(shared_elem2);
id_ac1(id_ac1==centerElement_id) = [];
id_ad1 = ac(shared_elem3);
id_ad1(id_ad1==centerElement_id) = [];

%%
ab1 = find(pairs(:,1)==a(1));
ab2 = find(pairs(:,2)==a(1));
ab = sort([ab1;ab2]);

ac1 = find(pairs(:,1)==a(2));
ac2 = find(pairs(:,2)==a(2));
ac = sort([ac1;ac2]);

ad1 = find(pairs(:,1)==a(3));
ad2 = find(pairs(:,2)==a(3));
ad = sort([ad1;ad2]);

shared_elem1 = ismember(ab,ac);
shared_elem2 = ismember(ab,ad);
shared_elem3 = ismember(ac,ad);

id_ab2 = ab(shared_elem1);
id_ac2 = ab(shared_elem2);
id_ad2 = ac(shared_elem3);

neighMatrixFrac = [id_ab2;id_ac2;id_ad2];

neighElement_id = [id_ab1;id_ac1;id_ad1];


% compare element
diff = ismember(tri_nodes(neighElement_id,:), pairs(neighMatrixFrac,:));
diff = double(diff);
diff = sum(diff,2);
diff = find(diff==2);

neighElement_id(diff) = [];
neighMatrixFracture_id = neighMatrixFrac + num_matrix;
%%
neighAll_id = [neighElement_id; neighMatrixFracture_id];