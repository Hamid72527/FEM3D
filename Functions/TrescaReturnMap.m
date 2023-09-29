function [eStress,eBStress,eSigmay,DGABAR,DDSDDE] = ...
           TrescaReturnMap(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel)
%------------------------------------------------------------------------
%  Purpose: 
%     STRESS UPDATE PROCEDURE FOR TRESCA TYPE ELASTO-PLASTIC MATERIAL WITH
%     ISOTROPIC/Kinematic HARDENING:
%     IMPLICIT ELASTIC PREDICTOR/RETURN MAPPING ALGORITHM.
%  Synopsis:
%     [eStress,eBStress,eSigmay,DGABAR,DDSDDE] = ...
%     TrescaReturnMap(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel)
%
%  Variable Description:
%     eStress - Stress vector
%     eBStress - Back stress vector
%     eSigmay - the yield stress
%     DGABAR - delta landa
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
ePEEQN=ePEEQ;
%**************************************************************************

eBStress=eBStressN;
DStTrial=Constants.DEV*StTrial;
PStTrial=(Constants.m')*StTrial/3;

eDBStress=Constants.DEV*eBStress;

% Spectral decomposition of the elastic trial deviatoric stress
[V,D] = eig(mkmatrixS(DStTrial));
PSTRS1=D(3,3);PSTRS2=D(2,2);PSTRS3=D(1,1);
e1=V(1:3,3);e2=V(1:3,2);e3=V(1:3,1);

% Spectral decomposition of the trial deviatoric back stress
[V,D] = eig(mkmatrixS(eDBStress));
PBSTRS1=D(3,3);PBSTRS2=D(2,2);PBSTRS3=D(1,1);
eB1=V(1:3,3);eB2=V(1:3,2);eB3=V(1:3,1);

SHMAXT=abs((PSTRS1-PBSTRS1)-(PSTRS3-PBSTRS3));
f_eq=SHMAXT-eSigmay;
 
% identify possible two-vector return: right or left of main plane
SCAPRD=(PSTRS1-PBSTRS1)+(PSTRS3-PBSTRS3)-2*(PSTRS2-PBSTRS2);
if SCAPRD>=0
    RIGHT=1; %.TRUE.
else
    RIGHT=0; %.FALSE.
end
%******************************************************
% Return mapping
%******************************************************
% Apply one-vector return mapping first (return to main plane)
Constants.TWOVEC=0; %.FALSE.
MXITER=500;
DGABAR=0;
for NRITER=1:MXITER
    % Compute residual derivative
    DENOM=-4*(MatProps.gmodule)-(MatProps.beta)*(MatProps.hmodule)-2*(2/3)*(MatProps.hmodule)*(1-(MatProps.beta));
    % Compute Newton-Raphson increment and update variable DGAMA
    dp=-f_eq/DENOM;
    DGABAR=DGABAR+dp;
    ePEEQ=ePEEQN+DGABAR;
    
    if isequal(MatProps.linhrd,"no")
        [eSigmay,MatProps.hmodule]=plfun(ePEEQ);
    else
        eSigmay=eSigmay+(MatProps.beta)*dp*(MatProps.hmodule);
    end
    
    % Compute new residual
    SHMAX=SHMAXT-4*(MatProps.gmodule)*DGABAR-(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*2*DGABAR;
    f_eq=SHMAX-eSigmay;
    
    % Check convergence
    RESNOR=abs(f_eq/eSigmay);
    if(RESNOR<=(1e-10))
        % Check validity of one-vector return
        S1=PSTRS1-2*(MatProps.gmodule)*DGABAR;
        S2=PSTRS2;
        S3=PSTRS3+2*(MatProps.gmodule)*DGABAR;
        BS1=PBSTRS1+(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*1*DGABAR;
        BS2=PBSTRS2;
        BS3=PBSTRS3+(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*-1*DGABAR;
        %DELTA=max([abs(S1),abs(S2),abs(S3)])*1.e-10;
        and((S1-BS1)>=(S2-BS2),(S2-BS2)>=(S3-BS3));
        %and((S1-BS1)+DELTA>=(S2-BS2),(S2-BS2)+DELTA>=(S3-BS3))
        if and((S1-BS1)>=(S2-BS2),(S2-BS2)>=(S3-BS3))
            % converged stress is in the same sextant as trial stress -> 1-vector
            % return is valid. Update principal deviatoric stresses
            PSTRS1=S1;
            PSTRS3=S3;
            PBSTRS1=BS1;
            PBSTRS3=BS3;
            break;
        else
            Constants.TWOVEC=1; %.TRUE.
            break
            % 1-vector return is not valid - go to two-vector procedure
        end
    end
    if NRITER==MXITER
        msg = 'failure of stress update procedure, ERROR-WE0001';
        error(msg)
    end
end
% ------------------------------------------------------------------
% Apply two-vector return mapping (return to corner - right or left)
% ------------------------------------------------------------------
if Constants.TWOVEC==1
    DGAMA=0;
    DGAMB=0;
    eSigmay=lastIncData.SigmayN(int*iel,1);
    SHMXTA=(PSTRS1-PBSTRS1)-(PSTRS3-PBSTRS3);
    if(RIGHT)
        SHMXTB=(PSTRS1-PBSTRS1)-(PSTRS2-PBSTRS2);
    else
        SHMXTB=(PSTRS2-PBSTRS2)-(PSTRS3-PBSTRS3);
    end
    PHIA=SHMXTA-eSigmay;
    PHIB=SHMXTB-eSigmay;
    % Start Newton-Raphson iterations
    for NRITER=1:MXITER
        % Compute residual derivative
        DRVAA=-(4*(MatProps.gmodule))-(MatProps.beta)*(MatProps.hmodule) -2*(2/3)*(MatProps.hmodule)*(1-(MatProps.beta));
        DRVAB=-(2*(MatProps.gmodule))-(MatProps.beta)*(MatProps.hmodule) -(2/3)*(MatProps.hmodule)*(1-(MatProps.beta));
        DRVBA=-(2*(MatProps.gmodule))-(MatProps.beta)*(MatProps.hmodule) -(2/3)*(MatProps.hmodule)*(1-(MatProps.beta));
        DRVBB=-(4*(MatProps.gmodule))-(MatProps.beta)*(MatProps.hmodule) -2*(2/3)*(MatProps.hmodule)*(1-(MatProps.beta));
        % Compute Newton-Raphson increment and update variables DGAMA and DGAMB
        R1DDET=1/(DRVAA*DRVBB-DRVAB*DRVBA);
        DDGAMA=(-DRVBB*PHIA+DRVAB*PHIB)*R1DDET;
        DDGAMB=(DRVBA*PHIA-DRVAA*PHIB)*R1DDET;
        DGAMA=DGAMA+DDGAMA;
        DGAMB=DGAMB+DDGAMB;
        % Compute new residual
        DGABAR=DGAMA+DGAMB;
        ePEEQ=ePEEQN+DGABAR;

        if isequal(MatProps.linhrd,"no")
            [eSigmay,MatProps.hmodule]=plfun(ePEEQ);
        else
            eSigmay=lastIncData.SigmayN(int*iel,1)+(MatProps.beta)*DGABAR*(MatProps.hmodule);
        end

        PHIA=SHMXTA-(2*(MatProps.gmodule))*(2*DGAMA+DGAMB)-(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*(2*DGAMA+DGAMB)-eSigmay;
        PHIB=SHMXTB-(2*(MatProps.gmodule))*(DGAMA+2*DGAMB)-(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*(DGAMA+2*DGAMB)-eSigmay;
        % Check convergence
        RESNOR=(abs(PHIA)+abs(PHIB))/eSigmay;
        if(RESNOR<=(1e-10))
            % Update EPBAR and principal deviatoric stresses
            if(RIGHT)
                PSTRS1=PSTRS1-(2*(MatProps.gmodule))*(DGAMA+DGAMB);
                PSTRS3=PSTRS3+(2*(MatProps.gmodule))*DGAMA;
                PSTRS2=PSTRS2+(2*(MatProps.gmodule))*DGAMB;
                % -------------------------------------
                PBSTRS1=PBSTRS1+(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*1*(DGAMA+DGAMB);
                PBSTRS2=PBSTRS2+(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*-1*DGAMB;
                PBSTRS3=PBSTRS3+(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*-1*DGAMA;
                break;
            else
                PSTRS1=PSTRS1-(2*(MatProps.gmodule))*DGAMA;
                PSTRS3=PSTRS3+(2*(MatProps.gmodule))*(DGAMA+DGAMB);
                PSTRS2=PSTRS2-(2*(MatProps.gmodule))*DGAMB;
                % -------------------------------------
                PBSTRS1=PBSTRS1+(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*1*DGAMA;
                PBSTRS2=PBSTRS2+(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*1*DGAMB;
                PBSTRS3=PBSTRS3+(1-(MatProps.beta))*(2/3)*(MatProps.hmodule)*-1*(DGAMA+DGAMB);
                break;
            end
        end
    end
    if NRITER==MXITER
        msg = 'failure of stress update procedure, ERROR-WE0002';
        error(msg)
    end
end

%------------------------------------
% update stress components
% ------------------------
eStress=(PSTRS1+PStTrial)*(e1*(e1'))+(PSTRS2+PStTrial)*(e2*(e2'))...
    +(PSTRS3+PStTrial)*(e3*(e3'));
eBStress=(PBSTRS1)*(eB1*(eB1'))+(PBSTRS2)*(eB2*(eB2'))...
    +(PBSTRS3)*(eB3*(eB3'));

eStress=mkvectorS(eStress);
eBStress=mkvectorS(eBStress);
%------------------------------------
if Constants.TWOVEC==0
    Constants.a=e1*e1'-e3*e3';
elseif RIGHT
    a2=e1*e1'-e2*e2';a1=e1*e1'-e3*e3';
    Constants.a=(DGAMA*a1+DGAMB*a2)/DGABAR;
else
    a2=e2*e2'-e3*e3';a1=e1*e1'-e3*e3';
    Constants.a=(DGAMA*a1+DGAMB*a2)/DGABAR;
end
Constants.a=mkvectorE(Constants.a);
DDSDDE=0.757*(Constants.matmtx-(1/((Constants.a')*Constants.matmtx*Constants.a))*Constants.matmtx*Constants.a*(Constants.a')*Constants.matmtx)+0.243*Constants.matmtx;

return
end