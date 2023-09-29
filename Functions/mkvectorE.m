function V=mkvectorE(M)
        V=[M(1,1);M(2,2);M(3,3);M(1,2)+M(2,1);M(2,3)+M(3,2);M(1,3)+M(3,1)];
    return
end