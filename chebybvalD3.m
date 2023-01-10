function[vector]=chebybvalD3(columns,eval)

vector=zeros(1,columns);
 
if eval==0

% build vector T = [ 1/2  0  -1  0  1  -1  0  1  ... ]:

T=zeros(1,columns);

    for i=1:4:columns
        T(i)=1;
    end

    for i=2:2:columns
        T(i)=0;
    end

    for i=3:4:columns
        T(i)=-1;
    end

T(1)=0.5;

% initialize auxiliary matrix:

M=zeros(columns,columns);
kk=0;

    for i=1:2:columns
        for j=3+kk:2:columns
            M(j,i)=T(i)*((j+1)^2-(i-1)^2)*((j-1)^2-(i-1)^2);
        end
        kk=kk+2;
    end

% build vector whose components are matrix M rows sums:

X=zeros(1,columns);

    for i=1:columns
        X(i)=sum(M(i,:));
    end

% build the actual bc vector for D^3W(0):


    for i=4:2:columns
        vector(i)=(i-1)/4*X(i-1);
    end

elseif eval==1

    for i=1:columns
        vector(i)=(i-1)^2*((i-1)^2-1)*((i-1)^2-4)/15;
    end
    
end

end
