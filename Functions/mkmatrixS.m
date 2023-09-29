function M=mkmatrixS(V)
    if or(length(V)==6,length(V)==9)
        M=[V(1) V(4) V(6);
           V(4) V(2) V(5);
           V(6) V(5) V(3)];
    else
        disp("Wrong Input!");
    end
    return
end