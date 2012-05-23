%Initialize the grid and the number of days
row = 20;
column = 20;
days = 365;

N = row*column;
colloc_data=zeros(N*days,3);
linearMatrix = zeros(row*column,3);
delta = 0.0001;

index=1;
for i=1:row
    for j=1:column
        linearMatrix(index,1) = i;
        linearMatrix(index,2) = j;
        linearMatrix(index,3) = index;
        index = index+1;
    end
end

sparseIndexX = zeros(N,1);
sparseIndexY = zeros(N,1);
sparseIndexValue = zeros(N,1);


index = 1;
for j=1:N
    x = linearMatrix(j,1);
    y = linearMatrix(j,2);
    for k=1:N
        x1 = linearMatrix(k,1);
        y1 = linearMatrix(k,2);

        if(j == k)
            %diagonal element
            sparseIndexX(index,1) = j;
            sparseIndexY(index,1) = k;


            % For the four corner elements
            if j == 1 || j == N || (x == 1 && y == column) || (x == row && y == 1)
                sparseIndexValue(index,1) = 8;
            elseif (x == 1 && y == 2)  || (x == 1 && y == column-1) || (x == 2 && y == 1)  || (x == 2 && y == column) || + ...
               (x == row-1 && y == 1) ||  (x == row-1 && y == column) || (x == row && y == 2) || (x == row && y == column-1)
                % For the elements at the edge but not the corner
                sparseIndexValue(index,1) = 11;
            elseif (x==2 && y==2) || (x==2 && y==column-1) || (x==row-1 && y==2) || (x==row-1 && y==column-1)
                sparseIndexValue(index,1) = 15;
            elseif (x==2 && y>2) || (x==2 && y<column-1) || (x==row-1 && y>2) || (x==row-1 && y<column-1) || (x>2 && y==2) || (x<row-1 && y==2) || (x>2 && y==column-1) || (x<row-1 && y==column-1)
                sparseIndexValue(index,1) = 19;
            elseif x > 2 && y > 2 && x < row-1 && y < column-1
                % For the elements in the middle
                sparseIndexValue(index,1) = 24;
            
            else
                sparseIndexValue(index,1) = 14;
            end
            
            index = index + 1;
        else
            %if the elements are adjacent 
            if (abs(x-x1) <= 2 && abs(y-y1) == 0) || (abs(x-x1) == 0 && abs(y-y1) <= 2) || (abs(x-x1) <= 2 && abs(y-y1) <= 2)
                sparseIndexX(index,1) = j;
                sparseIndexY(index,1) = k;
                sparseIndexValue(index,1) = -1;
                index = index + 1;
            end
        end
    end
end

%Generate the precision matrix for spatial correlation
QSpatial = sparse(sparseIndexX', sparseIndexY', sparseIndexValue', N, N);

beta=10;
Q=beta*QSpatial;

%Make the matrix strictly diagonally dominant so that it is positive definite and
%invertible
for i=1:N
    Q(i,i) = Q(i,i) + delta;
end

k=1;

%for k=1:days
    disp(k);
    v=zeros(N,1);
    
    %Get the cholesky decomposition
    %Step 1
    L=chol(Q,'lower');
    
    %Sample from the normal distribution
    %Step 2
    z=randn(1,N);
    z=z';
    
    %Back substitution method to generate the sample
    %Step 3
    for i=N:-1:1
        total = 0;
        if i~= N
            for j=(i+1):N
                total=total + L(j,i) * v(j,1);
            end
        end
        v(i,1) = 1/L(i,i) * (z(i,1) - total);
    end
    
    %Step 4
    x=v;
    
    start_index = (k-1)*N+1;
    end_index = (k-1)*N+N;
    colloc_data(start_index:end_index,3) = x;
    
%end
x=reshape(x,20,20);
imagesc(x);
colorbar;  