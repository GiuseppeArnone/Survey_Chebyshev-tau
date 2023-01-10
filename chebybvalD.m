function[vector]=chebybvalD(columns,eval)

vector=zeros(1,columns);
 
if eval==0

    for i=2:4:columns
        vector(i)=i-1;
    end

    for i=4:4:columns
        vector(i)=-i+1;
    end

elseif eval==1

    for i=1:columns
        vector(i)=(i-1)^2;
    end
    
end

end
