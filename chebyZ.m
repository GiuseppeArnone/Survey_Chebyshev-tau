function[matrixout]=chebyZ(rows,columns)

matrix=zeros(rows,columns);

for i=1:rows
    matrix(i,i+1)=0.5;
end

for i=1:columns
    matrix(i+1,i)=0.5;
end

    matrix(2,1)=1;

matrixout=matrix(1:rows,1:columns);

end



