function [a,hessian]=ahessCal(eStress,eBStress,eSigmay,PlasticPotential)
    DEV=(1/3)*[2 -1 -1 0 0 0;
          -1  2 -1 0 0 0;
          -1 -1  2 0 0 0;
           0  0  0 3 0 0;
           0  0  0 0 3 0;
           0  0  0 0 0 3];
    A=[-ones(3)+3*eye(3) zeros(3);zeros(3) 6*eye(3)];
    L=[eye(3) zeros(3);zeros(3) 2*eye(3)];
    %------------------------------------
    eDStress=DEV*eStress;
    eDBStress=DEV*eBStress;
    %------------------------------------
    if isequal(PlasticPotential,"Tresca")
        [V,~]= eig(mkmatrixS(eStress-eBStress));
        e1=V(1:3,3);e3=V(1:3,1);
        a=mkvectorE(e1*e1'-e3*e3');
        hessian=zeros(6,6);
    else
        eqStress=sqrt((3/2)*(eDStress-eDBStress)'*L*(eDStress-eDBStress));
        a=(3/(2*eSigmay))*L*(eDStress-eDBStress);
        hessian=A/(2*eqStress) - (1/eqStress)*a*(a');
    end
    return
end
