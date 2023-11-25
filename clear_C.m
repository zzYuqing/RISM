function C1=clear_C(C)
    C1=C;
    U_C=unique(C);
    num=length(U_C);
    for i=1:num
        points=find(C==U_C(i));
        C1(points)=i;
    end
end