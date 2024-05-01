function [row,col] = submultiple(N)
    row_init = ceil(sqrt(N));
    row_end = floor(row_init/2);
    row_vec = row_init:-1:row_end;
    col_vec = zeros(size(row_vec));
    res_vec = zeros(size(row_vec));
    for n = 1:length(row_vec)
        col_vec(n) = ceil(N/row_vec(n));
        res_vec(n) = row_vec(n)*col_vec(n)-N;
    end
    n = find(res_vec == min(res_vec),1);
    row = row_vec(n);
    col = col_vec(n);
end