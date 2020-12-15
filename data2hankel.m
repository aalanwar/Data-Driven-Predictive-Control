function H = data2hankel(data,num_block_rows)
%data =  (size of data entries) x (length of data stream)
%num_block_rows = number of block rows wanted in Hankel matrix
dim = size(data,1);
num_data_pts = size(data,2);
H = zeros(dim*(num_block_rows), num_data_pts-num_block_rows+1);

for i = 1:num_block_rows
    for j = 1:num_data_pts-num_block_rows+1
        H(dim*(i-1)+1:dim*i,j) = data(:,i+j-1);
    end
end

end