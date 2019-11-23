function [center_line, line_id, area1D, pairs] = centerLines(num_triangles, tri_nodes, pos_nodes, num_id)
% Finding the center of the lines and its area

pairs = zeros(3*num_triangles,2);
line = max(num_id);
id = 1;
for i=1:num_triangles
    for j=1:2
        k = j + 1;
        while k < 4
            pairs(id,:) = [tri_nodes(i,j), tri_nodes(i,k)];
            k = k + 1;
            id = id + 1;
        end
    end
end

o = 1;
while o < length(pairs)
    rem_same_pairs = pairs(o,:);
    
    ab(1) = rem_same_pairs(1);
    ab(2) = rem_same_pairs(2);
    finding_sim_val1 = find((pairs(:,1)==ab(1)) .* (pairs(:,2)==ab(2)));
    
    ab(1) = rem_same_pairs(2);
    ab(2) = rem_same_pairs(1);
    finding_sim_val2 = find((pairs(:,1)==ab(1)) .* (pairs(:,2)==ab(2)));
    pairs([finding_sim_val1(2:end),finding_sim_val2],:) = [];
    o = o + 1;
end

center_line = zeros(length(pairs),2);
line_id = zeros(length(pairs),1);
area1D = zeros(length(pairs),1);
for i=1:length(center_line)
    center_line(i,1) = mean(pos_nodes(pairs(i,:),1));
    center_line(i,2) = mean(pos_nodes(pairs(i,:),2));
    l1 = pos_nodes(pairs(i,1),:);
    l2 = pos_nodes(pairs(i,2),:);
    area1D(i,1) = pdist([l1;l2],'euclidean');
    line_id(i,1) = line + i;
end


