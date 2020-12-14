function [ U ] = arrange(U)
%arrange the u vector recursively in order to incoporate the nested 
%structure
U_mat = U;
[nrow,ncol] = size(U_mat);
if min(nrow,ncol)<=2
    U = reshape(U_mat',[],1);
    return
else
    cutting_row = round(ceil(nrow/2)+1-mod(nrow,2));
    cutting_col = round(ceil(ncol/2)+1-mod(ncol,2));
    %defind subblock of the matrix U_mat
    U11 = U_mat(1:cutting_row-1,1:cutting_col-1);
    U21 = U_mat(cutting_row+1:nrow,1:cutting_col-1);
    U12 = U_mat(1:cutting_row-1,cutting_col+1:ncol);
    U22 = U_mat(cutting_row+1:nrow,cutting_col+1:ncol);
    U = [arrange(U11);arrange(U21);U_mat(cutting_row,1:cutting_col-1)';...
        arrange(U12);arrange(U22);U_mat(cutting_row,cutting_col+1:ncol)';...
        U_mat(:,cutting_col)];
end

