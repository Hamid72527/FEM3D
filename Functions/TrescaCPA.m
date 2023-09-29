function [eStress,eBStress,eStranP,eSigmay,delp,DDSDDE] = ...
           TrescaCPA(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel)
%------------------------------------------------------------------------
%  Purpose: 
%     STRESS UPDATE PROCEDURE FOR TRESCA TYPE ELASTO-PLASTIC MATERIAL WITH
%     ISOTROPIC/Kinematic HARDENING:
%     EXPLICIT ELASTIC PREDICTOR/Cutting Plane ALGORITHM.
%  Synopsis:
%     [eStress,eBStress,eSigmay,DGABAR,DDSDDE] = ...
%     TrescaCPA(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel)
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
[Constants.a,Constants.hessian]=ahessCal(StTrial,eBStressN,eSigmay,MatProps.PlasticPotential);
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
Constants.TWOVEC=0;
while(1)
    i=i+1;
    eStress=eStress+deStress;            %1
    eBStress=eBStress+deBStress;         %2
    delp=delp+dp;                        %3
    eStranP=eStranP+deStranP;            %6
    ePEEQ=ePEEQ+dp;

    if isequal(MatProps.linhrd,"no")
        [eSigmay,MatProps.hmodule]=plfun(ePEEQ);
    else
        eSigmay=eSigmay+(MatProps.beta)*dp*(MatProps.hmodule);     %5
    end
    %------------------------------------    
    D = eig(mkmatrixS(eStress-eBStress));
    S1=D(3);S2=D(2);S3=D(1);
    eStress_eq=abs(S1-S3);
    %-----------------------------------
    %     check for yeild
    %-----------------------------------
    f_eq=eStress_eq-eSigmay;
    if norm(f_eq)<=0.001
        cond=and(S1>=S2,S2>=S3);
        if cond~=1
            Constants.TWOVEC=1;
            msg = '1-vector return is not valid , ERROR-WE0001';
            error(msg)
        end
        break;
    else
        if MatProps.kinhrd=="Prager"
            dp = f_eq/(Constants.a'*Constants.matmtx*Constants.a...
                +Constants.a'*(1-MatProps.beta)*(2/3)*MatProps.hmodule*Constants.invL*Constants.a+MatProps.beta*MatProps.hmodule);
            deBStress=(1-MatProps.beta)*(2/3)*MatProps.hmodule*Constants.invL*Constants.a*dp;
        elseif MatProps.kinhrd=="Ziegler"
            dp = f_eq/(Constants.a'*Constants.matmtx*Constants.a...
                +Constants.a'*(1-MatProps.beta)*MatProps.hmodule/eSigmay*(eStress-eBStress)+MatProps.beta*MatProps.hmodule);
            deBStress=(1-MatProps.beta)*MatProps.hmodule*(eStress-eBStress)*dp/eSigmay;
        end
        deStress=-Constants.matmtx*Constants.a*dp;
        deStranP=Constants.a*dp;

    end
end
%------------------------------------
%DDSDDE=Constants.matmtx;
%DDSDDE=Constants.matmtx-(1/((Constants.a')*Constants.matmtx*(Constants.a) + Constants.hmodule))*Constants.matmtx*(Constants.a)*(Constants.a')*Constants.matmtx;
DDSDDE=0.757*(Constants.matmtx-(1/((Constants.a')*Constants.matmtx*Constants.a))...
    *Constants.matmtx*Constants.a*(Constants.a')*Constants.matmtx)+0.243*Constants.matmtx;
%**************************************************************************
return
end