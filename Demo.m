%The small demo example
clear
for N = [128]
    x = 0:1/(N+1):1;
    y = 0:1/(N+1):1;
    [X,Y] = meshgrid(x,y);
    f_mat = sin(pi*X).*sin(pi*Y); 
    f_mat = f_mat(2:N+1,2:N+1);
    %Matrix A
    Z = zeros(N,N);
    Z((N+1:N+1:N^2-1)) = ones(N-1,1);
    A = -2*eye(N)+Z+Z';
    clear Z;
    A = (kron(eye(N),A)+kron(A,eye(N)));
    
    index = arrange(reshape((1:N^2),[],N)');
    A = (A(index,index));
    [~,invindex] = sort(index);
    
    [L,U] = LU_ND(A,0);
    f_mat = reshape(f_mat,[],1);
    
    
    u_vec = U\(L\f_mat(index))./(N+1)^2;
    
    surf(reshape(u_vec,[],N));
    %surf(reshape(u_vec(invindex),[],N))
end