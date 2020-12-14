%Comparing Multifrontal solver with Dense linear algebra
clear; 

for N = [5]
%N = 21;
    Z = zeros(N,N);
    Z((N+1:N+1:N^2-1)) = ones(N-1,1);
    A = -2*eye(N)+Z+Z';
    clear Z;
    A = (kron(eye(N),A)+kron(A,eye(N)));
    %imagesc((1:N^2),(1:N^2),A);
    spy(A)

    %Solving transformation index
    index = arrange(reshape((1:N^2),[],N)');
    %Modity matrix A accordingly
    A = sparse(A(index,index));
    %imagesc((1:N^2),(1:N^2),A);
   % spy(A)

    tic
    [L,U] = LU_ND(A,0);
    
    toc
    fprintf('The error of LU(multifrontal server) when N = %d is %d \n',...
        N,norm(A-L*U,inf)) 
    
end
%[L,U] = lu(A);

%-------------------Dense linear algebra-----------------------------------
for N = [5]%,10,15,20,25,30,35,45,55,65,75,85,95,115,128]
    Z = zeros(N,N);
    Z((N+1:N+1:N^2-1)) = ones(N-1,1);
    A = -2*eye(N)+Z+Z';
    clear Z;
    A = kron(eye(N),A)+kron(A,eye(N));
    %imagesc((1:N^2),(1:N^2),A);
    %spy(A)

    %Solving transformation index
    index = arrange(reshape((1:N^2),[],N)');
    %Modity matrix A accordingly
    A = A(index,index);
    %imagesc((1:N^2),(1:N^2),A);
    %spy(A)

    tic
    [L,U] = lu(A);
    
    toc
    fprintf('The error of LU(Dense linear algebra) when N = %d is %d \n',...
        N,norm(A-L*U,inf)) 
end

