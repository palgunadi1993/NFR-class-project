function [neigh_fractureID, neigh_matrixID] = fracture_matrix_fracture(ID_fracture, line_ID, pairs, tri_nodes)
% Find neighbor of fracture: matrices and fractures

ifmf = find(line_ID==ID_fracture);
a = pairs(ifmf,:);

% Fracture-matrix
ab1 = find(tri_nodes(:,1)==a(1));
ab2 = find(tri_nodes(:,2)==a(1));
ab3 = find(tri_nodes(:,3)==a(1));
ab = sort([ab1;ab2;ab3]);

ac1 = find(tri_nodes(:,1)==a(2));
ac2 = find(tri_nodes(:,2)==a(2));
ac3 = find(tri_nodes(:,3)==a(2));
ac = sort([ac1;ac2;ac3]);

shared_elem1 = ismember(ab,ac);

neigh_matrixID = ab(shared_elem1);

% fracture-fracture
fa1 = find(pairs(:,1)==a(1));
fa2 = find(pairs(:,2)==a(1));
fa = sort([fa1;fa2]);

fb1 = find(pairs(:,1)==a(2));
fb2 = find(pairs(:,2)==a(2));
fb = sort([fb1;fb2]);

fa(fa==ifmf) = [];
fb(fb==ifmf) = [];

neigh_fractureID = [fa;fb] + length(tri_nodes);