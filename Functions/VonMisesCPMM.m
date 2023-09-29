function [eStress,eBStress,eStranP,eSigmay,delp,DDSDDE] = ...
           VonMisesCPMM(lastIncData,MatProps,Constants,DSTRAN,StTrial,ePEEQ,int,iel)       
%------------------------------------------------------------------------
%  Purpose: 
%     STRESS UPDATE PROCEDURE FOR VON-MISES TYPE ELASTO-PLASTIC MATERIAL WITH
%     ISOTROPIC/Kinematic HARDENING:
%     EXPLICIT ELASTIC PREDICTOR/CLOSEST POINT ALGORITHM.
%  Synopsis:
%     [eStress,eBStress,eSigmay,DGABAR,DDSDDE] = ...
%     VonMisesCPMM(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel)
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
eStressN=lastIncData.StressN(int*iel,:)';
eBStressN=lastIncData.BStressN(int*iel,:)';
%eStranN=lastIncData.StranN(int*iel,:)';
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
    eBStress=eBStress+deBStress;         %2
    delp=delp+dp;                        %3
    delS=eStress-(eStressN);             %4
    delBS=eBStress-(eBStressN);          %5
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
    [Constants.a,Constants.hessian]=ahessCal(eStress,eBStress,eSigmay,MatProps.PlasticPotential);
    %------------------------------------    
    R_k=[Constants.matmtx\delS-DSTRAN+Constants.a*delp;
         Constants.hrdmtx\delBS-Constants.a*delp*(1-MatProps.beta)];
    f_eq=eStress_eq-eSigmay;
    %-----------------------------------
    %     check for yeild
    %-----------------------------------
    if i==1
        R_tmp=R_k;
    end
    cond1=or(norm(R_k(1:6))<=0.001*norm(R_tmp(1:6)),norm(R_k(1:6))<=0.001);
    cond2=or(norm(R_k(7:12))<=0.001*norm(R_tmp(7:12)),norm(R_k(7:12))<=0.001);
    cond3=and(cond1,cond2);
    if and(norm(f_eq)<=0.001,cond3)
        break;
    else
        Am=[inv(Constants.matmtx) zeros(6);zeros(6) inv(Constants.hrdmtx)]...
        +delp*[Constants.hessian, -Constants.hessian*(1-MatProps.beta);
              -Constants.hessian, Constants.hessian*(1-MatProps.beta)];
        invAm=eye(12)/Am;

        H_star=[Constants.a;-Constants.a]'*invAm*[Constants.a;-Constants.a]+MatProps.beta*MatProps.hmodule;
        dp=(-1/H_star)*([Constants.a;-Constants.a]'*invAm*R_k-f_eq);
        delk=-invAm*R_k+(1/H_star)*invAm*[Constants.a;-Constants.a]...
            *([Constants.a;-Constants.a]'*invAm*R_k-f_eq);
        deStress=delk(1:6);
        deBStress=delk(7:12);
        deStranP=dp*Constants.a;
    end
end
%------------------------------------
K_tmp=invAm-(1/H_star)*invAm*[Constants.a;-Constants.a]*[Constants.a;-Constants.a]'*invAm;
DDSDDE=(Constants.I_star')*K_tmp*(Constants.I_star);
return
end