function interpol = divided_difference_q5(T, Y, t_1)

% construct divided difference table
j = length(T);

table = zeros(j, j); % create a table of 0s with size jxj

table(:,1) = Y; % fill the first column of the table with the values of Y

% use for loop to iterate through the table
for z = 2:j 
    for i = z:j 
        table(i,z) = (table(i,z-1) - table(i-1,z-1)) / (T(i) - T(i+1-z));
    end
end 

% evaluate divided difference
m = size(table);

interpol = zeros(size(t_1)); % initialise sum of interpolating polynomial
for k = 1:m % for loop to iterate through table(k, k)
    P = ones(size(t_1)); % initialise product for interpol
    for i = 1:k-1 % for loop to iterate through T
        P = P .* (t_1-T(i)); % multiply the next factor
    end
    interpol = interpol + table(k,k) * P; % add the next term
end

end