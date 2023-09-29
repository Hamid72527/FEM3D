function [matmtrx]=fematiso(iopt,elastic,poisson)

%------------------------------------------------------------------------
%  Purpose:
%     determine the constitutive equation for isotropic material
%
%  Synopsis:
%     [matmtrx]=fematiso(iopt,elastic,poisson) 
%
%  Variable Description:
%     elastic - elastic modulus
%     poisson - Poisson's ratio   
%     iopt=1 - plane stress analysis
%     iopt=2 - plane strain analysis
%     iopt=3 - axisymmetric analysis
%     iopt=4 - three dimensional analysis - vector form
%     iopt=5 - three dimensional analysis - index form
%------------------------------------------------------------------------

 if iopt==1        % plane stress
   matmtrx= elastic/(1-poisson*poisson)* ...
   [1  poisson 0; ...
   poisson  1  0; ...
   0  0  (1-poisson)/2];

 elseif iopt==2     % plane strain
   matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson 0; 
   poisson  (1-poisson)  0;
   0  0  (1-2*poisson)/2];

 elseif iopt==3     % axisymmetry
   matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson  0; 
   poisson  (1-poisson)   poisson  0;
   poisson  poisson  (1-poisson)   0;
   0    0    0   (1-2*poisson)/2];
 
 elseif iopt==4    % three-dimension - vector form(6*6)
   matmtrx= elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2];

 elseif iopt==5    % three-dimension - index form(3*3*3*3)
   mu=elastic/(2*(1+poisson));
   landa=-(2/3)*mu+elastic/(3*(1-2*poisson));

   delta=diag(ones(3,1));
   matmtrx=zeros(3,3,3,3);

   for i=1:3
       for j=1:3
           for k=1:3
               for l=1:3
                   matmtrx(i,j,k,l)=landa*delta(i,j)*delta(k,l)+mu*...
                       (delta(i,k)*delta(j,l)+delta(i,l)*delta(j,k));
               end
           end
       end
   end
 elseif iopt==6 % three-dimension - vector form(9*9)
   P=0.5*[2 0 0 0 0 0;
       0 2 0 0 0 0;
       0 0 2 0 0 0;
       0 0 0 1 0 0;
       0 0 0 1 0 0;
       0 0 0 0 1 0;
       0 0 0 0 1 0;
       0 0 0 0 0 1;
       0 0 0 0 0 1];  
   matmtrx= P*(elastic/((1+poisson)*(1-2*poisson))* ...
   [(1-poisson)  poisson  poisson   0   0    0; 
   poisson  (1-poisson)   poisson   0   0    0;
   poisson  poisson  (1-poisson)    0   0    0;
   0    0    0    (1-2*poisson)/2   0    0;
   0    0    0    0    (1-2*poisson)/2   0;
   0    0    0    0    0   (1-2*poisson)/2])*P';
 end