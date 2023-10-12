function [p_B_7, p_B_6, p_B_5, p_B_4, p_B_3, p_B_2, p_B_1, p_B_0] = p_B_committor(cell_N1, cell_N2, cell_N3,cell_N4,cell_N5,cell_N6, cell_N7, k0, k1,k2,k3,k4,k5,k6, k7)

mat_N1 = cell2mat(cell_N1); %number of succesfull trajectories initiating from all points at lambda_1
mat_N2 = cell2mat(cell_N2);%number of succesfull trajectories initiating from all points at lambda_2
mat_N3 = cell2mat(cell_N3);%number of succesfull trajectories initiating from all points at lambda_3
mat_N4 = cell2mat(cell_N4);%number of succesfull trajectories initiating from all points at lambda_4
mat_N5 = cell2mat(cell_N5);%number of succesfull trajectories initiating from all points at lambda_5
mat_N6 = cell2mat(cell_N6); %number of succesfull trajectories initiating from all points at lambda_6
mat_N7 = cell2mat(cell_N7); %number of succesfull trajectories initiating from all points at lambda_6

p_B_6 = zeros(length(mat_N6),1);
p_B_5 = zeros(length(mat_N5),1);
p_B_4 = zeros(length(mat_N4),1);
p_B_3 = zeros(length(mat_N3),1);
p_B_2 = zeros(length(mat_N2),1);
p_B_1 = zeros(length(mat_N1),1);

%p_B_6 = mat_N6./k6;
p_B_7 = mat_N7./k7;

cumulative_mat_N6 = cumsum(mat_N6);
cumulative_mat_N5 = cumsum(mat_N5);
cumulative_mat_N4 = cumsum(mat_N4);
cumulative_mat_N3 = cumsum(mat_N3);
cumulative_mat_N2 = cumsum(mat_N2);
cumulative_mat_N1 = cumsum(mat_N1);

p_B_6(1) = sum(p_B_7(1:cumulative_mat_N6(1)))/k6;
p_B_5(1) = sum(p_B_6(1:cumulative_mat_N5(1)))/k5;
p_B_4(1) = sum(p_B_5(1:cumulative_mat_N4(1)))/k4;
p_B_3(1) = sum(p_B_4(1:cumulative_mat_N3(1)))/k3;
p_B_2(1) = sum(p_B_3(1:cumulative_mat_N2(1)))/k2;
p_B_1(1) = sum(p_B_2(1:cumulative_mat_N1(1)))/k1;


for j = 2:length(mat_N6)

    %Basically, we are developing a mapping scheme;
    % For the last crossing point at lambda_5, we know that 3 trjacetories
    % reached lambda_6. Hence, on the cumulative vector for lambda_5, for slicing from p_B_6,  the
    % upper idx will simply be the last element of cumsum vector for lambda_5, whereas, the lower
    % idx will be the penultiumate (i.e., j-1) element + 1 (e.g., last
    % index is 734469, and penultimate is 734466. Since 3 trajectories
    % reached lambda_6, from p_B_6, we need to slice 734467: 734469
    
    idx_lower = cumulative_mat_N6(j-1)+1; 
    idx_upper = cumulative_mat_N6(j);
    p_B_6(j) = sum(p_B_7(idx_lower:idx_upper))/k6;
 
end

for j = 2:length(mat_N5)

    %Basically, we are developing a mapping scheme;
    % For the last crossing point at lambda_5, we know that 3 trjacetories
    % reached lambda_6. Hence, on the cumulative vector for lambda_5, for slicing from p_B_6,  the
    % upper idx will simply be the last element of cumsum vector for lambda_5, whereas, the lower
    % idx will be the penultiumate (i.e., j-1) element + 1 (e.g., last
    % index is 734469, and penultimate is 734466. Since 3 trajectories
    % reached lambda_6, from p_B_6, we need to slice 734467: 734469
     
    idx_lower = cumulative_mat_N5(j-1)+1; 
    idx_upper = cumulative_mat_N5(j);
    p_B_5(j) = sum(p_B_6(idx_lower:idx_upper))/k5;
 
end


for j = 2:length(mat_N4)
     
    p_B_4(j) = sum(p_B_5(cumulative_mat_N4(j-1)+1:cumulative_mat_N4(j)))/k4;
    
end

for j = 2:length(mat_N3)
     
    p_B_3(j) = sum(p_B_4(cumulative_mat_N3(j-1)+1:cumulative_mat_N3(j)))/k3;
    
end

for j = 2:length(mat_N2)
     
    p_B_2(j) = sum(p_B_3(cumulative_mat_N2(j-1)+1:cumulative_mat_N2(j)))/k2;
    
end

for j = 2:length(mat_N1)
     
    p_B_1(j) = sum(p_B_2(cumulative_mat_N1(j-1)+1:cumulative_mat_N1(j)))/k1;
    
end

p_B_0 = sum(p_B_1)/k0;
    
end

%[p_B_7, p_B_6, p_B_5, p_B_4, p_B_3, p_B_2, p_B_1, p_B_0] = p_B_committor(cell_N1, cell_N2, cell_N3,cell_N4,cell_N5,cell_N6, cell_N7, k0, k1,k2,k3,k4,k5,k6, k7)
