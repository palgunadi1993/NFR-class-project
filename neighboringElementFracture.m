function [neighMatrixFracture_id] = neighboringElementFracture(tri_nodes,pairs, centerElement_id)
% Matrix-fracture

a = tri_nodes(centerElement_id,:);
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

id_ab = ab(shared_elem1);
id_ac = ab(shared_elem2);
id_ad = ac(shared_elem3);

neighMatrixFracture_id = [id_ab;id_ac;id_ad];