function [eStress,eBStress,eStranP,eSigmay,delp,DDSDDE] = ...
           TrescaCPMM(lastIncData,MatProps,Constants,DSTRAN,StTrial,ePEEQ,int,iel)      
%------------------------------------------------------------------------
%  Purpose: 
%     STRESS UPDATE PROCEDURE FOR TRESCA TYPE ELASTO-PLASTIC MATERIAL WITH
%     ISOTROPIC/Kinematic HARDENING:
%     IMPLICIT ELASTIC PREDICTOR/ClOSEST POINT ALGORITHM.
%  Synopsis:
%     [eStress,eBStress,eSigmay,DGABAR,DDSDDE] = ...
%     TrescaCPMM(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel)
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
eStranPN=lastIncData.StranPN(int*iel,:)';
%******************************************************
[Constants.a,Constants.hessian]=ahessCal(StTrial,eBStressN,eSigmay,MatProps.PlasticPotential);

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
    % Prager Kinematic Hardening rule
    Cp=2*MatProps.hmodule/3;
    Constants.hrdmtx=Cp*eye(9);
    Constants.hrdmtx=(Constants.P1')*Constants.hrdmtx*(Constants.P1);
    %------------------------------------
    D2=eig(mkmatrixS(eBStress));
    D1=eig(mkmatrixS(eStress));
    D=D1-D2;
    S1=D(3);S2=D(2);S3=D(1);
    eStress_eq=abs(S1-S3);
    
%     D = eig(mkmatrixS(eStress-eBStress));
%     S1=D(3);S2=D(2);S3=D(1);
%     eStress_eq=abs(S1-S3);
    %------------------------------------
    R_k=[Constants.matmtx\delS-DSTRAN+Constants.a*delp;
         Constants.hrdmtx\delBS-Constants.a*delp*(1-MatProps.beta)];
    if i==1
        R_tmp=R_k;
    end
    f_eq=eStress_eq-eSigmay;
    %-----------------------------------
    %     check for yeild
    %-----------------------------------
    cond1=or(norm(R_k(1:6))<=0.001*norm(R_tmp(1:6)),norm(R_k(1:6))<=0.001);
    cond2=or(norm(R_k(7:12))<=0.001*norm(R_tmp(7:12)),norm(R_k(7:12))<=0.001);
    cond3=and(cond1,cond2);
    if and(and(norm(f_eq)<=0.001,cond3),i~=1)
        cond=and(S1>=S2,S2>=S3);
        if cond~=1
            Constants.TWOVEC=1;
            if (S1-S2-eSigmay>0)
                Constants.RIGHT=1;
            else
                Constants.RIGHT=0;
            end
        end
        break;
    else
        Am=[inv(Constants.matmtx) zeros(6);zeros(6) inv(Constants.hrdmtx)]...
        +delp*[Constants.hessian, -Constants.hessian*(1-MatProps.beta);
              -Constants.hessian, Constants.hessian*(1-MatProps.beta)];
        invAm=eye(12)/Am;

        H_star=[Constants.a;-Constants.a]'*invAm*[Constants.a;-Constants.a]+(MatProps.beta)*(MatProps.hmodule);
        dp=(-1/H_star)*([Constants.a;-Constants.a]'*invAm*R_k-f_eq);
        delk=-invAm*R_k+(1/H_star)*invAm*[Constants.a;-Constants.a]...
            *([Constants.a;-Constants.a]'*invAm*R_k-f_eq);
        deStress=delk(1:6);
        deBStress=delk(7:12);
        deStranP=dp*Constants.a;
    end
end
K_tmp=invAm-(1/H_star)*invAm*[Constants.a;-Constants.a]*[Constants.a;-Constants.a]'*invAm;
% ------------------------------------------------------------------
% Apply two-vector return mapping (return to corner - right or left)
% ------------------------------------------------------------------
if Constants.TWOVEC==1
    if MatProps.beta~=0
        msg='! This method only can solve problems with kinematic hardening materials !, ERROR-WE0002';
        error(msg);
    end
    %**********************************************************************
    % Extract the last force increment's(n-1) results
    eSigmay=lastIncData.SigmayN(int*iel,1);
    eStressN=lastIncData.StressN(int*iel,:)';
    eBStressN=lastIncData.BStressN(int*iel,:)';
    eStranPN=lastIncData.StranPN(int*iel,:)';
    %**********************************************************************
    eStranP=eStranPN;    %1
    delp=zeros(2,1);     %2
    eBStress=eBStressN;  %3
    eStress=StTrial;     %4
    %5
    deStress=zeros(6,1);
    deBStress=zeros(6,1);
    dp=zeros(2,1);
    deStranP=zeros(6,1);
    %**********************************************************************
    Constants.hessian=zeros(6,6);
    [V,~]= eig(mkmatrixS(eStress-eBStress));
    e1=V(1:3,3);e2=V(1:3,2);e3=V(1:3,1);
    if Constants.RIGHT
        a1=mkvectorE(e1*e1'-e3*e3');
        a2=mkvectorE(e1*e1'-e2*e2');
    else
        a1=mkvectorE(e1*e1'-e3*e3');
        a2=mkvectorE(e2*e2'-e3*e3');
    end
    Constants.a=[a1,a2];
    %**********************************************************************
    i=0;
    while(1)
        i=i+1;
        eStress=eStress+deStress;            %1
        eBStress=eBStress+deBStress;         %2
        delp=delp+dp;                        %3
        delS=eStress-(eStressN);             %4
        delBS=eBStress-(eBStressN);          %5
        eStranP=eStranP+deStranP;            %6
        ePEEQ=ePEEQ+dp(1)+dp(2);
        
        if isequal(MatProps.linhrd,"no")
            MatProps.hmodule=-5e10*ePEEQ+16e9;
        end

        eSigmay=eSigmay+(MatProps.beta)*(dp(1)+dp(2))*(MatProps.hmodule);     %5
        %------------------------------------
        % Prager Kinematic Hardening rule
        Cp=2*(MatProps.hmodule)/3;
        Constants.hrdmtx=Cp*eye(9);
        Constants.hrdmtx=(Constants.P1')*Constants.hrdmtx*(Constants.P1);
        %------------------------------------
        D2=eig(mkmatrixS(eBStress));
        D1=eig(mkmatrixS(eStress));
        D=D1-D2;
        S1=D(3);S2=D(2);S3=D(1);
        eStress_eq1=abs(S1-S3);
        if Constants.RIGHT
            eStress_eq2=abs(S1-S2);
        else
            eStress_eq2=abs(S2-S3);
        end
        f_eq1=eStress_eq1-eSigmay;
        f_eq2=eStress_eq2-eSigmay;
        f_eq=[f_eq1;f_eq2];
        
        %------------------------------------
        R_k=[Constants.matmtx\delS-DSTRAN+(Constants.a)*delp;
             Constants.hrdmtx\delBS-(Constants.a)*delp*(1-MatProps.beta)];
        if i==1
            R_tmp=R_k;
        end
        %-----------------------------------
        %     check for yeild
        %-----------------------------------
        cond1=or(norm(R_k(1:6))<=0.001*norm(R_tmp(1:6)),norm(R_k(1:6))<=0.001);
        cond2=or(norm(R_k(7:12))<=0.001*norm(R_tmp(7:12)),norm(R_k(7:12))<=0.001);
        cond3=and(cond1,cond2);
        cond4=and(norm(f_eq1)<=0.001,norm(f_eq2)<=0.001);
        if and(and(cond4,cond3),i~=1)
            break;
        else
            Am=[inv(Constants.matmtx) zeros(6);zeros(6) inv(Constants.hrdmtx)]...
            +delp(1)*[Constants.hessian, -Constants.hessian*(1-MatProps.beta);
                  -Constants.hessian, Constants.hessian*(1-MatProps.beta)]...
            +delp(2)*[Constants.hessian, -Constants.hessian*(1-MatProps.beta);
                  -Constants.hessian, Constants.hessian*(1-MatProps.beta)];  
            invAm=eye(12)/Am;
            G=[Constants.a;-Constants.a]'*invAm*[Constants.a;-Constants.a];
            dp=(-inv(G))*([Constants.a;-Constants.a]'*invAm*R_k-f_eq);
            delk=-invAm*R_k+invAm*[Constants.a;-Constants.a]*(eye(2)/G)...
                *([Constants.a;-Constants.a]'*invAm*R_k-f_eq);
            deStress=delk(1:6);
            deBStress=delk(7:12);
            deStranP=Constants.a*dp;
        end
    end
    delp=delp(1)+delp(2);
    K_tmp=invAm-invAm*[Constants.a;-Constants.a]*(eye(2)/G)*[Constants.a;-Constants.a]'*invAm;
end
%------------------------------------
DDSDDE=(Constants.I_star')*K_tmp*(Constants.I_star);
%**************************************************************************
return
end