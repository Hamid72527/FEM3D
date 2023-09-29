function [eStress,eBStress,eStranP,eSigmay,delp,DDSDDE] = ...
           VonMisesCPA(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel)      
%------------------------------------------------------------------------
%  Purpose: 
%     STRESS UPDATE PROCEDURE FOR VON-MISES TYPE ELASTO-PLASTIC MATERIAL WITH
%     ISOTROPIC/Kinematic HARDENING:
%     EXPLICIT ELASTIC PREDICTOR/CUTTING PLANE ALGORITHM.
%  Synopsis:
%     [eStress,eBStress,eSigmay,DGABAR,DDSDDE] = ...
%     VonMisesCPA(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel)
%
%  Variable Description:
%     eStress - Stress vector
%     eBStress - Back stress vector
%     eSigmay - the yield stress
%     delp - delta landa
%     DDSDDE - Tangent modulus
%     lastIncData - an array contains the last increments datas
%     MatProps - an array contains the material properties
%     Constants - an array contains some constant matrices
%     StTrial - the trial stress vector
%     ePEEQ - the accumulated plastic strain
%     int - the gauss points number which is been under analysis
%     iel - the element number which is been under analysis
%------------------------------------------------------------------------ 
       
%**************************************************************************
% Extract the last force increment's(n-1) results
eSigmay=lastIncData.SigmayN(int*iel,1);
eBStressN=lastIncData.BStressN(int*iel,:)';
eStranPN=lastIncData.StranPN(int*iel,:)';
%******************************************************

eStranP=eStranPN;    %1
delp=zeros(1,1);     %2
eBStress=eBStressN;  %3
eStress=StTrial;     %4
%5
deStress=zeros(6,1);
deBStress=zeros(6,1);
dp=zeros(1,1);
deStranP=zeros(6,1);
%------------------------------------
i=0;
while(1)
    i=i+1;
    eStress=eStress+deStress;            %1
    eDStress=Constants.DEV*eStress;
    eBStress=eBStress+deBStress;         %2
    eDBStress=Constants.DEV*eBStress;
    delp=delp+dp;                        %3
    eStranP=eStranP+deStranP;            %6
    ePEEQ=ePEEQ+dp;

    if isequal(MatProps.linhrd,"no")
        [eSigmay,MatProps.hmodule]=plfun(ePEEQ);
    else
        eSigmay=eSigmay+(MatProps.beta)*dp*(MatProps.hmodule);     %5
    end
    %------------------------------------
    if MatProps.kinhrd=="Prager"
        Cp=2*MatProps.hmodule/3;
        Constants.hrdmtx=Cp*eye(9);
    elseif MatProps.kinhrd=="Ziegler"
        Cz=MatProps.hmodule/eSigmay;
        Constants.hrdmtx=(Cz*2*eSigmay/3)*eye(9);
    end
    Constants.hrdmtx=(Constants.P1')*Constants.hrdmtx*(Constants.P1);
    %------------------------------------
    % calculate the Effective Trial Stress(eStress_eq)
    [eStress_eq]=eqStressCal(eStress,eBStress,MatProps.PlasticPotential);
    Constants.a=(3/2/eStress_eq)*Constants.L*(eDStress-eDBStress);
    %------------------------------------    
    f_eq=eStress_eq-eSigmay;
    %-----------------------------------
    %     check for yeild
    %-----------------------------------
    if norm(f_eq)<=0.001
        break;
    else
        dp=f_eq/((Constants.a'*Constants.matmtx*Constants.a)...
            +(-(1-MatProps.beta)*Constants.a'*(Constants.hrdmtx)*Constants.a)+MatProps.beta*MatProps.hmodule);
        deBStress=(1-MatProps.beta)*dp*(Constants.hrdmtx)*(Constants.a);
        deStress=-Constants.matmtx*Constants.a*dp;
        deStranP=Constants.a*dp;
    end
end
%------------------------------------
DDSDDE=0.757*(Constants.matmtx-(1/((Constants.a')*Constants.matmtx*Constants.a))...
    *Constants.matmtx*Constants.a*(Constants.a')*Constants.matmtx)+0.243*Constants.matmtx;
return
end