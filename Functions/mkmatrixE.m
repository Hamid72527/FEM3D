function M=mkmatrixE(V)
    if or(length(V)==6,length(V)==9)
        M=[V(1)   ,V(4)/2 ,V(6)/2;
           V(4)/2 ,V(2)   ,V(5)/2;
           V(6)/2 ,V(5)/2 ,V(3)];
    else
        disp("Wrong Input!");
    end
    return
end