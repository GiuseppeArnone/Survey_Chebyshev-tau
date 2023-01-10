function[vector]=chebybval(columns,eval)

vector=zeros(1,columns);
 
if eval==0

    for i=1:4:columns
        vector(i)=1;
    end

    for i=2:2:columns
        vector(i)=0;
    end

    for i=3:4:columns
        vector(i)=-1;
    end

elseif eval==1

    for i=1:columns
        vector(i)=1;
    end
    
end

end
