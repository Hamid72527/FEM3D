%##########################################################################                                                          
%##---------------------- Computational Plasticity code -----------------##
%##----------------------------- Hamid Alijai ---------------------------##
%########################################################################## 

clc
clear
close all
feval('setpath')
format
tic
%--------------------------------------------------------------------------
Model  = createpde('structural','static-solid');

%OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
%  input data for control parameters
%OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

RadialReturn="Tresca-CPA"; %"Tresca-ReturnMap" "Tresca-CPMM" "Tresca-CPA"
                           %"VonMises-CPA" "VonMises-CPMM"
%------------------------------------
eltype="linear"; %'quadratic' or 'linear'
seedsize=3;
a_coeff=3.0e+05; 

%------------------------------------
MatProps.emodule=2e10;            % elastic modulus [Kg/m^3]
MatProps.hmodule=2/3*1e10;        % plastic tangent modulus [Kg/m^2](define if hardening law in linear)
MatProps.poisson=0.3;             % Poisson's ratio
MatProps.kmodule=MatProps.emodule/(3*(1-2*MatProps.poisson)); % bulk modulus [Kg/m^3]
MatProps.gmodule=MatProps.emodule/(2*(1+MatProps.poisson));   % shear modulus (G) [Kg/m^3]
MatProps.Sigma0 = 2000e4;  % Yield Stress [Kg/m^2]
MatProps.beta=0.5;           % Beta (an indicator to determine the type of Hardening rule)
MatProps.kinhrd="Prager";  % "Prager" or "Ziegler"
MatProps.linhrd="yes"; % "yes" or "no"
%------------------------------------
dt=0.001;
t_start=0.0;
t_end=1.0;
incnum=round((t_end-t_start)/0.01);
%------------------------------------
nelplot=1;
%------------------------------------
% filename="Part.stl";
% importGeometry(Model,filename);
%------------------------------------
W=1 ;D=1 ; H=4;
Model.Geometry = multicuboid(W,H,D);
%OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO
%OOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOOO

%==========================================================================
%  Mesh
%====================================
mesh=generateMesh(Model,'Hmax',seedsize,'GeometricOrder',eltype);
% pdegplot(Model,'FaceLabels','on','FaceAlpha',0.5)
% figure
% pdeplot3D(mesh,'NodeLabels','on');

%==========================================================================
%  Set another parameters
%====================================
if eltype=="linear"
    nnel=4;                  % number of nodes per element
    ngl=1;                   % number of GL points per element
elseif eltype=="quadratic"
    nnel=10;                  % number of nodes per element
    ngl=4;                    % number of GL points per element
end
ndof=3;                  % number of dofs per node    
edof=nnel*ndof;          % degrees of freedom per element
NDI=3;  % Number of direct stress components at this point.
NSHR=3; % Number of engineering shear stress components at this point.
NTENS=NDI+NSHR; % Size of the stress or strain component array (NDI + NSHR).
%------------------------------------                           
if or(isequal(RadialReturn,"VonMises-CPA"),isequal(RadialReturn,"VonMises-CPMM"))
   MatProps.PlasticPotential="VonMises"; 
else
   MatProps.PlasticPotential="Tresca";
end
%==========================================================================
%  input data for nodal coordinate values
%  gcoord(i,j) where (i-1)->node no. and j->x or y or z
%====================================
gcoord=Model.Mesh.Nodes';
tmp=size(gcoord);
nnode=tmp(1);            % total number of nodes in system
sdof=nnode*ndof;         % total system dofs  

%==========================================================================
%  input data for nodal connectivity for each element
%  nodes(i,j) where i-> element no. and j-> connected nodes
%====================================
nodes=Model.Mesh.Elements';
tmp=size(nodes);
nel=tmp(1);                   % number of elements

%==========================================================================
%  input data for boundary conditions
%   bcdof -> a vector containing dofs associated with boundary conditions     
%   bcval -> a vector containing boundary condition values associated with    
%           the dofs in 'bcdof'                                              
%====================================
% figure
% pdegplot(Model,'FaceLabels','on','FaceAlpha',0.5)
NBC = findNodes(mesh,'region','Face',6);
Nf=findNodes(mesh,'region','Face',4);
% hold on
% p1=plot3(mesh.Nodes(1,NBC),mesh.Nodes(2,NBC),mesh.Nodes(3,NBC),'ok','MarkerFaceColor','b');  
% p2=plot3(mesh.Nodes(1,Nf),mesh.Nodes(2,Nf),mesh.Nodes(3,Nf),'ok','MarkerFaceColor','r');  
% legend([p1 p2],"BC-Fixed","Load");
%-----------------------------------------
bcdof=[NBC*3,NBC*3-1,NBC*3-2]; %X,Y and Z direction of Face ID 6 have been constrained
% for tetrahedron

% pdeplot3D(mesh,'NodeLabels','on');
% Nf=3;NBC=[1 2 4];
% bcdof=[1 2 3 4 5 6 10 11 12]; % nine dofs are constrained
%-----------------------------------------
bcval=zeros(1,length(bcdof));   % whose described values are 0 

%==========================================================================
%  initialization of matrices and vectors
%====================================
index=zeros(edof,1);    % a vector containing system dofs associated with each element
kinmtx=zeros(6,edof);   % kinematic matrix (B)
xcoord=zeros(1,4);      % vector containing x coordination of element's nodes
ycoord=zeros(1,4);      % vector containing y coordination of element's nodes
zcoord=zeros(1,4);      % vector containing z coordination of element's nodes
nd=zeros(1,4);          % vector containing node numbers of an element
edelu=zeros(edof,1);    % element nodal differential displacement vector
eStress=zeros(NTENS,1); % the stress vector at an integration point
eBStress=zeros(NTENS,1); % the back stress vector at an integration point
%-----------------------------------------
ff1=zeros(sdof,1);  % external system force vector 
ff2=zeros(sdof,1);  % internal system force vector 
kk=zeros(sdof,sdof);   % system matrix
DDSDDE=zeros(NTENS,NTENS); % Jacobian matrix of the constitutive Model (C_ep tangent)
displ=zeros(sdof,1);       % vector containing the nodal displacements
Stress=zeros(ngl*nel,NTENS);  % matrix containing stress components
BStress=zeros(ngl*nel,NTENS); % matrix containing back stress components
Strain=zeros(ngl*nel,NTENS);  % matrix containing strain components
StranP=zeros(ngl*nel,NTENS);  % matrix containing strain components
Sigmay=MatProps.Sigma0*ones(ngl*nel,1); % vector containing the yield stress at integration points
PEEQ=zeros(nel,ngl);
PEEQ_tmp=zeros(nel,ngl);
%-----------------------------------------
lastIncData.SigmayN=Sigmay;                % the yield stress vector at last load increment
lastIncData.StressN=zeros(ngl*nel,NTENS);  % the stress vector at last load increment
lastIncData.BStressN=zeros(ngl*nel,NTENS); % the back stress vector at last load increment
lastIncData.StranN=zeros(ngl*nel,NTENS);  % the strain vector at last load increment
lastIncData.StranPN=zeros(ngl*nel,NTENS);  % the plastic strain vector at last load increment

Results.Stress.Tresca=zeros(nel,ngl,incnum);
Results.Stress.TrescaNodal=zeros(nnode,1);

Results.Stress.Pressure=zeros(nel,ngl,incnum);
Results.Stress.PressureNodal=zeros(nnode,1);

Results.Stress.Mises=zeros(nel,ngl,incnum);
Results.Stress.MisesNodal=zeros(nnode,1);

Results.Strain.PEEQ=zeros(nel,ngl,incnum);
Results.Strain.PEEQNodal=zeros(nnode,1);

Results.Strain.MaxPrincipal=zeros(nel,ngl,incnum);
Results.Strain.MaxPrincipalNodal=zeros(nnode,1);

Results.Strain.MidPrincipal=zeros(nel,ngl,incnum);
Results.Strain.MinPrincipal=zeros(nel,ngl,incnum);
Results.Displacement.Magnitude=zeros(nnode,1);
Results.Displacement.ux=zeros(nnode,1);
Results.Displacement.uy=zeros(nnode,1);
Results.Displacement.uz=zeros(nnode,1);
%-----------------------------------------
% Define some matrices
Constants.P1=0.5*[2 0 0 0 0 0;
       0 2 0 0 0 0;
       0 0 2 0 0 0;
       0 0 0 1 0 0;
       0 0 0 1 0 0;
       0 0 0 0 1 0;
       0 0 0 0 1 0;
       0 0 0 0 0 1;
       0 0 0 0 0 1];
   
Constants.P2=[1 0 0 0 0 0;
    0 1 0 0 0 0;
    0 0 1 0 0 0;
    0 0 0 1 0 0;
    0 0 0 1 0 0;
    0 0 0 0 1 0;
    0 0 0 0 1 0;
    0 0 0 0 0 1;
    0 0 0 0 0 1];

Constants.L=[eye(3) zeros(3);zeros(3) 2*eye(3)];
Constants.invL=eye(6)/Constants.L;

Constants.A=[-ones(3)+3*eye(3) zeros(3);zeros(3) 6*eye(3)];
Constants.I_star=[eye(6);zeros(6)];

Constants.M2=(1/3)* [2 -1 -1 0 0 0;
          -1  2 -1 0 0 0;
          -1 -1  2 0 0 0;
           0  0  0 6 0 0;
           0  0  0 0 6 0;
           0  0  0 0 0 6];

% define Deviatoric matirx(DEV)
% dstress=DEV*stress, dstrain=DEV*strain
Constants.DEV=(1/3)*[2 -1 -1 0 0 0;
          -1  2 -1 0 0 0;
          -1 -1  2 0 0 0;
           0  0  0 3 0 0;
           0  0  0 0 3 0;
           0  0  0 0 0 3];
Constants.m=[1;1;1;0;0;0];

%-----------------------------------------
if Model.Mesh.GeometricOrder =="quadratic"
    alpha=0.30872448322;
    beta=1.928411633;
    gama=0.8098434003;
    NodesCoords=[-alpha -alpha -alpha;
                  beta  -alpha -alpha;
                 -alpha  beta  -alpha;
                 -alpha -alpha  beta;
                  gama  -alpha -alpha;
                  gama   gama  -alpha;
                 -alpha  gama  -alpha;
                 -alpha -alpha  gama;
                  gama  -alpha  gama;
                 -alpha  gama   gama];
    Constants.RecoveryMatrix=zeros(10,4);
    for ii=1:10
        rvalue=NodesCoords(ii,1);svalue=NodesCoords(ii,2);tvalue=NodesCoords(ii,3);
        [shapet4,~,~,~]=feisot4(rvalue,svalue,tvalue);
        Constants.RecoveryMatrix(ii,:)=shapet4;
    end
else
    Constants.RecoveryMatrix=[1;1;1;1];
end

%==========================================================================
%  computation of element matrices and vectors and their assembly
%====================================

%   point3 -> matrix containing sampling points
%   weight3 -> matrix containing weighting coefficients
[point3,weight3]=feglqd3t(ngl);  % sampling points & weights

%-----------------------------------------
Constants.matmtx=fematiso(4,MatProps.emodule,MatProps.poisson);        % compute constitutive matrix
n=0;
l=0;
for t=t_start:dt:t_end
    n=n+1;
    delu=zeros(sdof,1);
    du=zeros(sdof,1);
    
    %-----------------------------------------
    %  force vector
    ff1=zeros(sdof,1);  % external system force vector
    %ff1(Nf*3-2)=a_coeff*forceX(t)/length(Nf);    % force applied at node 4 in x-axis
    %ff1(Nf*3-1)=a_coeff*forceY(t)/length(Nf);  % force applied at node 4 in y-axis
    ff1(Nf*3)=-a_coeff*forceZ(t)/1;         % force applied at node 4 in z-axis
    %**********************************************************************
    j=0;
    convergence1 = 0;
    while convergence1 == 0     % start LOOP newton-raphson
        j=j+1;
        delu=delu+du;

        ff2=zeros(sdof,1);     % external system force vector
        kk=zeros(sdof,sdof);   % system matrix

        for iel=1:nel           % loop for the total number of elements
            for i=1:nnel
                nd(i)=nodes(iel,i);       % extract connected node for (iel)-th element
                xcoord(i)=gcoord(nd(i),1);  % extract x value of the node
                ycoord(i)=gcoord(nd(i),2);  % extract y value of the node
                zcoord(i)=gcoord(nd(i),3);  % extract z value of the node
            end
            
            k=zeros(edof,edof);         % initialization of element matrix to zero
            f=zeros(edof,1);
            index=feeldof(nd,nnel,ndof);% extract system dofs associated with element
            
            for i=1:edof
                edelu(i,1)=delu(index(i));
            end
            %--------------------------------
            %  numerical integration
            %--------------------------------
            for int=1:ngl
                rvalue=point3(int,2);            % sampling point in r-axis
                svalue=point3(int,3);            % sampling point in s-axis
                tvalue=point3(int,4);            % sampling point in t-axis
                wtr=weight3(int);                % weight
                if eltype=="linear"
                    [shape,dhdr,dhds,dhdt]=feisot4(rvalue,svalue,tvalue);
                    % compute shape functions and derivatives at sampling point
                elseif eltype=="quadratic"
                    [shape,dhdr,dhds,dhdt]=feisot10(rvalue,svalue,tvalue);
                    % compute shape functions and derivatives at sampling point
                end
                jacob3=fejacob3(nnel,dhdr,dhds,dhdt,xcoord,ycoord,zcoord);  % compute Jacobian
                detjacob=det(jacob3);                 % determinant of Jacobian
                invjacob=inv(jacob3);                 % inverse of Jacobian matrix
                [dhdx,dhdy,dhdz]=federiv3(nnel,dhdr,dhds,dhdt,invjacob); % derivatives w.r.t. physical coordinate
                kinmtx=fekine3d(nnel,dhdx,dhdy,dhdz);          % compute kinematic matrix
                
                %----------------------------------------------------------
                % Element de-Assemble
                eStress=lastIncData.StressN(int*iel,:)';
                eBStress=lastIncData.BStressN(int*iel,:)';
                eDBStress=Constants.DEV*eBStress;
                edelu;
                eSigmay=lastIncData.SigmayN(int*iel,1);

                ePEEQ=PEEQ(iel,int);
                
                DSTRAN=kinmtx*edelu;
                Strain(int*iel,:)=lastIncData.StranN(int*iel,:)+DSTRAN';

                %----------------------------------------------------------
                % calculate the Trial Stress(STrial)
                StTrial=eStress+Constants.matmtx*DSTRAN;
 
                % calculate the Deviatoric Trial Strian(DStTrial)
                DStTrial=Constants.DEV*StTrial;
                
                % calculate the Effective Trial Stress(ESTrial) , Hessian and a 
                [ESTrial]=eqStressCal(StTrial,eBStress,MatProps.PlasticPotential);
                %-----------------------------------    
                %     check for yeild
                %-----------------------------------
                f_eq=ESTrial-eSigmay;
                if f_eq<=0.001
                    % update stress and tangent modulus
                    eStress=StTrial;
                    Stress(int*iel,:)=eStress';
                    BStress(int*iel,:)=eBStress';
                    DDSDDE=Constants.matmtx;
                else
                    switch RadialReturn
                        
                        %-Tresca-------------------------------------------
                        case "Tresca-CPA"
                            [eStress,eBStress,eStranP,eSigmay,delp,DDSDDE] = ...
                        TrescaCPA(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel);
                        StranP(int*iel,:) = eStranP';
                        case "Tresca-ReturnMap"
                            [eStress,eBStress,eSigmay,delp,DDSDDE] = ...
                        TrescaReturnMap(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel);
                        case "Tresca-CPMM"
                            [eStress,eBStress,eStranP,eSigmay,delp,DDSDDE] = ...
                        TrescaCPMM(lastIncData,MatProps,Constants,DSTRAN,StTrial,ePEEQ,int,iel);
                        StranP(int*iel,:) = eStranP';
                    
                        %- Von Mises --------------------------------------
                        case "VonMises-CPA"
                            [eStress,eBStress,eStranP,eSigmay,delp,DDSDDE] = ...
                        VonMisesCPA(lastIncData,MatProps,Constants,StTrial,ePEEQ,int,iel);
                        StranP(int*iel,:) = eStranP';
                        case "VonMises-CPMM"
                            [eStress,eBStress,eStranP,eSigmay,delp,DDSDDE] = ...
                        VonMisesCPMM(lastIncData,MatProps,Constants,DSTRAN,StTrial,ePEEQ,int,iel);
                        StranP(int*iel,:) = eStranP';
                        otherwise
                            msg="Error-Wrong Input for RadialReturn";
                            error(msg);
                    end
                    %------------------------------------------------------
                    Stress(int*iel,:)= eStress';
                    BStress(int*iel,:)= eBStress';
                    Sigmay(int*iel,1)=eSigmay;
                    PEEQ_tmp(iel,int)=delp+PEEQ(iel,int);
                end
                k=k+kinmtx'*DDSDDE*kinmtx*wtr*detjacob*1/6;    % element matrix                    
                f=f+kinmtx'*eStress*wtr*detjacob*1/6;
            end % end of GL loop
            kk=feasmbl1(kk,k,index);  % assemble element matrices
            ff2=feasmbl2(ff2,f,index); % assemble element vectors
        end % end of Element loop

        ff=ff1-ff2;
        [kk,ff]=feaplyc2(kk,ff,bcdof,bcval);
        du=kk\ff;
         % compute some variables and Check the convergance criteria
        Residual=norm(ff);
        cond1=or(norm(du)<=0.001*norm(delu),norm(du)<=0.001*norm(displ));
        if j<=2
            Residual_tmp=Residual;
        end
        cond2=or(Residual<=0.1*Residual_tmp,Residual<=0.1);
         if and(cond2,cond1)
            clc
            disp("================================");
            disp("Processing ... ");
            disp("================================");
            disp("Residual");
            disp(norm(ff));
            disp("norm(du) [nm]:");
            disp(norm(du)*1e9);
            disp("t is [ms]");
            disp(t*1000);
            disp("dt is [ms]:")
            disp(dt*1000);
            disp("iteration number is:");
            disp(j);
            disp("================================");
            
            convergence1 = 1;
            lastIncData.StressN=Stress;
            lastIncData.BStressN=BStress; 
            lastIncData.StranN=Strain;
            lastIncData.StranPN=StranP;
            lastIncData.SigmayN=Sigmay;
            PEEQ=PEEQ_tmp;
        end
    end % End of N-R Iteration [Output=delu]
    %**********************************************************************            
    displ=displ+delu;
    % Store everythings !!!
    if and(t~=0,mod(t,0.01)==0)
        l=l+1;
        for iel=1:nel
            for int=1:ngl
                eStrain=Strain(int*iel,:)';
                PSTRAN = eig(mkmatrixE(eStrain));
                Results.Strain.MaxPrincipal(iel,int,l)=PSTRAN(3);
                Results.Strain.MidPrincipal(iel,int,l)=PSTRAN(2);
                Results.Strain.MinPrincipal(iel,int,l)=PSTRAN(1);
                Results.Strain.PEEQ(iel,int,l)=PEEQ(iel,int);
                %----------------------------------------------------------                
                eStress=Stress(int*iel,:)';
                eBStress=BStress(int*iel,:)';
                [~,J2,J3]=stressinv(eStress);
                theta=acos(3*sqrt(3)/2*J3/(J2^1.5))/3;
                Results.Stress.Tresca(iel,int,l)=2*sqrt(J2)*sin(theta+pi/3);
                Results.Stress.Mises(iel,int,l)=sqrt(3*J2);
                Results.Stress.Pressure(iel,int,l)=(Constants.m')*eStress/3;
                %----------------------------------------------------------
            end
        end
    end
end

% -------------------------------------------------------------------------
close all
%pdeplot3D(Model.Mesh,'ElementLabels','on');
% tmp=size(Results.Stress.Mises);
% X=zeros(tmp(3),ngl);
% Y=zeros(tmp(3),ngl);
% for i=1:tmp(3)
%     for j=1:ngl
%         X(i,j)=Results.Strain.MaxPrincipal(nelplot,j,i);
%         Y(i,j)=Results.Stress.Mises(nelplot,j,i);
%     end
% end
% figure
% hold on
% for i=1:ngl
%     Xtmp=X(:,i);
%     Ytmp=Y(:,i);
%     plot(Xtmp,Ytmp)
% end
% xlabel("\bf \epsilon_{\rm eff}")
% ylabel("\bf \sigma_{\rm eff [Kg/m^2]}")
% title({'\sigma - \epsilon'; ['(for element number #',num2str(nelplot),')']},'Color','b')
% -------------------------------------------------------------------------

for i=1:(length(displ)/3)
    Results.Displacement.ux(i)=displ(3*i-2);
    Results.Displacement.uy(i)=displ(3*i-1);
    Results.Displacement.uz(i)=displ(3*i);
    Results.Displacement.Magnitude(i)=sqrt((displ(3*i)^2)+...
       (displ(3*i-1)^2)+(displ(3*i-2)^2));
end
% -------------------------------------------------------------------------
clc
disp("********************************************************************");
disp("----------------------------- Results ------------------------------");
disp("********************************************************************");
disp("Equivalent Tresca Stresses at GPs are [kg/m^2]:");
disp(Results.Stress.Tresca(:,:,end));
disp("Equivalent VonMises Stresses at GPs are [kg/m^2]:");
disp(Results.Stress.Mises(:,:,end));
disp("********************************************************************");
disp(" The Max displacements is [cm]: ");
disp(max(Results.Displacement.Magnitude)*100);
disp("********************************************************************");
toc
% -------------------------------------------------------------------------
% Interpolate Stresses to nodes
eTrescaNodal=zeros(nnel,1);
for iel=1:nel
    nd=nodes(iel,:);       % extract connected node for (iel)-th element
    % ----------------------------------
    eVMNodal=(Constants.RecoveryMatrix)*Results.Stress.Mises(iel,:,end)';
    eTrescaNodal=(Constants.RecoveryMatrix)*Results.Stress.Tresca(iel,:,end)';
    ePEEQ=(Constants.RecoveryMatrix)*Results.Strain.PEEQ(iel,:,end)';
    eStran=(Constants.RecoveryMatrix)*Results.Strain.MaxPrincipal(iel,:,end)';
    ePressure=(Constants.RecoveryMatrix)*Results.Stress.Pressure(iel,:,end)';
    % ----------------------------------
    for ii=1:nnel
        if (Results.Stress.MisesNodal(nd(ii))==0)
            Results.Stress.MisesNodal(nd(ii))=eVMNodal(ii);
        else
            Results.Stress.MisesNodal(nd(ii))=(eVMNodal(ii)...
                +Results.Stress.MisesNodal(nd(ii)))/2;
        end
        % ----------------------------------
        if (Results.Stress.TrescaNodal(nd(ii))==0)
            Results.Stress.TrescaNodal(nd(ii))=eTrescaNodal(ii);
        else
            Results.Stress.TrescaNodal(nd(ii))=(eTrescaNodal(ii)...
                +Results.Stress.TrescaNodal(nd(ii)))/2;
        end
        % ----------------------------------
        if (Results.Strain.PEEQNodal(nd(ii))==0)
            Results.Strain.PEEQNodal(nd(ii))=ePEEQ(ii);
        else
            Results.Strain.PEEQNodal(nd(ii))=(ePEEQ(ii)...
                +Results.Strain.PEEQNodal(nd(ii)))/2;
        end
        % ----------------------------------
        if (Results.Strain.MaxPrincipalNodal(nd(ii))==0)
            Results.Strain.MaxPrincipalNodal(nd(ii))=eStran(ii);
        else
            Results.Strain.MaxPrincipalNodal(nd(ii))=(eStran(ii)...
                +Results.Strain.MaxPrincipalNodal(nd(ii)))/2;
        end
        % ----------------------------------
        if (Results.Stress.PressureNodal(nd(ii))==0)
            Results.Stress.PressureNodal(nd(ii))=ePressure(ii);
        else
            Results.Stress.PressureNodal(nd(ii))=(ePressure(ii)...
                +Results.Stress.PressureNodal(nd(ii)))/2;
        end
    end
end
% -------------------------------------------------------------------------
figure('Name','Measured Data','Position',[300 50 750 600]);
subplot(2,2,1)
pdeplot3D(Model,'ColorMapData',Results.Stress.TrescaNodal,...
                'Deformation',Results.Displacement,'DeformationScaleFactor',25);
c = colorbar;
title('Tresca Stress')
c.Label.String = 'Tresca Stress (Kg/m^{2})';

subplot(2,2,2)
pdeplot3D(Model,'ColorMapData',Results.Stress.MisesNodal,...
                'Deformation',Results.Displacement,...
                'DeformationScaleFactor',25);
title('Mises Stress in the Beam for the Last Time-Step')
c = colorbar;
title('Mises Stress')
c.Label.String = 'Mises Stress (Kg/m^{2})';
% -------------------------------------------------------------------------
subplot(2,2,[3,4])
pdeplot3D(Model,'ColorMapData',Results.Displacement.Magnitude,...
                'Deformation',Results.Displacement,...
                'DeformationScaleFactor',25);
title('Displacement in the Beam for the Last Time-Step')
hold on
p1=plot3(mesh.Nodes(1,NBC),mesh.Nodes(2,NBC),mesh.Nodes(3,NBC),'ok','MarkerFaceColor','b');
legend(p1,"BC-Fixed",'Location','northwest');
c = colorbar;
c.Label.String = 'Displacement (m)';
% -------------------------------------------------------------------------



