function[matrixout]=chebyder(rows,columns)

matrix=zeros(rows,columns);

ck = 2;
for i=1:rows
    for j=i+1:2:columns
        matrix(i,j)=2*(j-1)/ck;
    end
ck=1;
end

matrixout=matrix;

end

