function [eqStress]=eqStressCal(eStress,eBStress,PlasticPotential)
%------------------------------------------------------------------------
%  Purpose:
%     to determine the equal stress in the framework of the Tresca criteria
%     or VonMises.   
%
%  Synopsis:
%     [eqStress]=eqStressCal(eStress,eBStress,PlasticPotential)  
%
%  Variable Description:
%     eStress - Stress vector
%     eBStress - Back stress vector
%     PlasticPotential -the yield criteria which is used, it could be
%     "Tresca" or "VonMises".
%------------------------------------------------------------------------    

    DEV=(1/3)*[2 -1 -1 0 0 0;
          -1  2 -1 0 0 0;
          -1 -1  2 0 0 0;
           0  0  0 3 0 0;
           0  0  0 0 3 0;
           0  0  0 0 0 3];
    L=[eye(3) zeros(3);zeros(3) 2*eye(3)];
    %------------------------------------
    eDStress=DEV*eStress;
    eDBStress=DEV*eBStress;
    %------------------------------------
    if isequal(PlasticPotential,"Tresca")
        PSTRS1= eig(mkmatrixS(eStress));
        PSTRS2= eig(mkmatrixS(eBStress));
        PSTRS= PSTRS1-PSTRS2;
        eqStress=abs(PSTRS(3)-PSTRS(1));
    else
        eqStress=sqrt((3/2)*(eDStress-eDBStress)'*L*(eDStress-eDBStress));
    end
    return
end
