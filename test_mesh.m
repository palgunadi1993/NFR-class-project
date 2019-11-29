clear 
close all
clc

load complex_mesh.mat

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

% an = [0.684079629	5.907990315;
% 1.127606192	6.004842615;
% 1.580544082	6.501210654;
% 1.960963264	6.86440678;
% 2.423246743	7.687651332;
% 2.985033909	8.462469734;
% 3.22966266	8.777239709;
% 3.736910149	9.285714286;
% 4.198864944	9.745762712];

% an = [2.976586723	9.128329298;
%     3.22966266	8.777239709;
% 3.374185137	8.474576271;
% 3.663339651	7.99031477;
% 4.060861372	7.251815981;
% 4.440338326	6.573849879;
% 4.747537607	6.02905569;
% 4.865009368	5.83535109;
% 5.126992648	5.326876513;
% 5.506502471	4.685230024;
% 5.804695803	4.188861985;
% 6.184391879	3.753026634;
% 6.500783364	3.365617433;
% 6.817317279	3.13559322;
% 7.124768551	2.869249395;
% 7.441258642	2.590799031;
% 7.866170719	2.118644068;
% 8.227800116	1.719128329;
% 8.489991564	1.440677966;
% 8.706967011	1.198547215];

% an = [1.349303737	0.98062954;
% 1.747932028	1.464891041;
% 2.110262619	1.840193705;
% 2.436339334	2.15496368;
% 2.744294589	2.445520581;
% 3.052359405	2.857142857;
% 3.351385404	3.280871671;
% 3.668478082	3.668280872;
% 3.867934658	4.06779661;
% 4.13985516	4.539951574;
% 4.493366056	5.169491525;
% 4.701905273	5.605326877;
% 4.865009368	5.83535109;
% 5.145913907	6.234866828;
% 5.725647234	6.840193705;
% 6.332497014	7.409200969;
% 6.930264153	7.94188862;
% 7.636716225	8.571428571;
% 7.953710298	8.849878935];

% an = [0.728890252	5.423728814;
% 1.271635643	5.157384988;
% 1.877674668	4.830508475;
% 2.302783956	4.576271186;
% 2.836567221	4.406779661;
% 3.487931809	4.16464891;
% 3.867934658	4.06779661;
% 4.211639806	3.861985472;
% 4.980476154	3.426150121;
% 5.695101509	3.08716707;
% 6.210631841	2.748184019;
% 6.518203631	2.615012107;
% 6.979457233	2.300242131;
% 7.540214521	1.937046005;
% 7.92007494	1.682808717];

% an = [2.272435441	1.041162228;
% 2.743133238	1.162227603;
% 3.268228282	1.392251816;
% 3.820417867	1.561743341;
% 4.318243073	1.658595642;
% 4.915648658	1.791767554;
% 5.594721331	2.167070218;
% 6.228468441	2.457627119;
% 6.518203631	2.615012107;
% 6.835055274	2.736077482;
% 7.124768551	2.869249395;
% 7.459774523	3.050847458;
% 7.939653567	3.317191283;
% 8.365179188	3.523002421;
% 8.645864604	3.680387409;
% 8.998926298	3.813559322];
% 
% an_u = [an(:,1)-0.0025 an(:,2) + 0.0025];
% an_d = [an(:,1)+0.0025 an(:,2) - 0.0025];
% 
% [in,on] = inpolygon(center_fracture(:,1),center_fracture(:,2), [an_u(:,1);flip(an_d(:,1))],[an_u(:,2);flip(an_d(:,2))]);
% inon = in | on; 
% idx = find(inon(:));
% xcoord = center_fracture(idx,1);
% ycoord = center_fracture(idx,2);
% fracture1 = fracture_id(idx,1);
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
% plot(an(:,1),an(:,2),'-g')
% plot(an_u(:,1),an_u(:,2),'-r')
% plot(an_d(:,1),an_d(:,2),'-r')
% plot(xcoord,ycoord, 'gp')
% plot([an_u(:,1);an_d(:,1)], [an_u(:,2);an_d(:,2)], '-r')
%%

% Find fracture
% x_bor = 5.0;
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
% y_bor = 5.0;
% x_bor = [2.0,8.0];
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

% x_bor = 5.0;
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
% c = a;

% Initiation of fracture
% frac = [161, 182, 201, 223, 251, 278, 299, 314];
% frac = [210, 231,250,273,306,325,344,188,214,244, 272,304,341,362,380];
% frac = [1148,1100,1050,1006,974,948,931,924,900,916,952,971,993,1025,1064,1119,1164,1215,1263,...
%     1282,1315,1321,1310,1300,1156,1083,1040,982,929,876,835,796,764,741];

% frac = [411,444,490,526,559];
% frac = c;% cmg

frac = [19369
19589
19832
20046
20287
20592
20851
21143
21432
21715
21992
22321
22630
22906
23189
23496
23795
24088
24437
24768
25143
25554
25942
26316
26697
27100
27472
27885
28322
28697
29101
29488
29851
30228
30594
30941
31282
31640
31982
32345
32754
33072
33459
33756
34036
34330
34573
34878
35129
35351
35582
35787
35984
36158
36327
36474
36643
27298
27304
27305
27320
27330
27344
27345
27388
27410
27416
27425
27429
27445
27463
27469
27498
27507
27513
27516
27517
27521
27553
27558
27574
27575
27577
27584
27585
27592
27605
27623
27624
27636
27673
27676
27689
27691
27698
27715
27717
27719
27724
27735
27743
27760
27762
27799
27812
27918
28031
28128
28249
28308
28419
28528
28627
28694
28811
28876
28970
29112
29169
29340
29476
29621
29710
29789
29872
29920
30027
30100
30213
30296
30395
30484
30580
30673
30772
30866
30931
31063
31138
31248
31309
31439
31502
31625
31742
31872
31919
31984
32128
32236
32319
32401
32511
32570
32686
32782
32845
32967
33084
33160
33238
33288
33299
33399
33428
15543
15578
15615
15655
15688
15739
15796
15846
15895
15969
16059
16118
16204
16286
16365
16460
16558
16657
16742
16844
16949
17054
17169
17295
17407
17531
17702
17855
18032
18168
18375
18576
18777
19047
19262
19524
19844
20119
20388
20692
20951
21191
21497
21732
21976
22298
22564
22818
23146
23436
23725
24040
24336
24754
25004
25440
25833
26150
26538
26900
27283
27692
28094
28548
28885
29303
29631
30053
30407
30743
31109
31461
31814
32147
32482
32883
33168
33497
33840
34108
34426
34707
34989
35242
35475
35679
35876
36078
36241
36388
36575
36721
36845
36965
37053
37168
37253
37349
37437
37504
37586
37657
37713
37777
37840
37907
37961
38021
38061
38106
38159
38197
38238
18440
18540
18603
18699
18751
18833
18885
18982
19036
19090
19194
19289
19342
19448
19519
19602
19651
19737
19809
19883
20006
20147
20290
20384
20548
20662
20811
20902
21034
21182
21256
21390
21525
21701
21869
22006
22134
22197
22251
22358
22402
22475
22608
22659
22735
22860
22949
23009
23110
23229
23299
23449
23522
23655
23771
23916
24015
24131
24160
24264
24333
24396
24472
24539
24640
24760
24938
25023
25158
25216
25241
25302
25373
25464
25488
25588
25674
25775
25870
25879
25982
26096
26121
26236
26274
26419
15879
15939
16023
16086
16148
16243
16324
16417
16511
16600
16699
16792
16894
16992
17120
17223
17353
17465
17584
17694
17839
17947
18079
18251
18423
18608
18813
18988
19187
19406
19664
19937
20213
20471
20778
21059
21329
21606
21882
22166
22454
22757
23091
23400
23714
24057
24342
24628
24962
25334
25623
25908
26188
26605
26927
27242
27550
28002
28345
28692
29063
29444
29811
30142
30479
30814
31127
31456
31785
32122
32456
32814
33106
33387
33641
33870
34146
34386
34678
34949];

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