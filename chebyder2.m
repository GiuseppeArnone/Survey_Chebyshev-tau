function[matrixout]=chebyder2(rows,columns)

matrix=zeros(rows,columns);

ck = 2;
for i=1:rows
    for j=i+2:2:columns
        matrix(i,j)=(j-1)*((j-1)^2-(i-1)^2)/ck;
    end
ck=1;
end

matrixout=matrix;

end



