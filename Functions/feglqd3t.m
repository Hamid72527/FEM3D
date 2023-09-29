function [point3,weight3]=feglqd3t(ngl)

%-------------------------------------------------------------------
%  Purpose:
%     determine the integration points and weighting coefficients of 
%     Gauss-Legendre quadrature for 3D tetrahedra integration
%
%  Synopsis:
%     [point3,weight3]=feglqd3t(ngl)
%
%  Variable Description:
%     ngl - number of integration points
%     point3 - vector containing integration points   
%     weight3 - vector containing weighting coefficients
%
%     Reference : https://www.sharcnet.ca/Software/Ansys/17.0/
%                 en-us/help/ans_thry/thy_et1.html
%-------------------------------------------------------------------

%  initialization

   point3=zeros(ngl,4);
   weight3=zeros(ngl,4);

 if ngl==1 % 1-point rule
     % Center Point	
     L1=.250000000000000;
     W1=1.000000000000000;
     point3=[L1,L1,L1,L1];
     weight3=W1;
 elseif ngl==4 % 4-point rule
     % Corner Points	
     L1=.585410196624968;
     L2=.138196601125010;
     W1=0.250000000000000;
     point3=[L1,L2,L2,L2
             L2,L1,L2,L2;
             L2,L2,L1,L2;
             L2,L2,L2,L1];
     weight3=[W1;W1;W1;W1];
 elseif ngl==5 % 5-point rule
     % Center Point
     L1=.250000000000000;
     W1=-0.800000000000000;
     % Corner Points	
     L2=.500000000000000;
     L3=.166666666666666;
     W2=0.450000000000000;
     point3=[L1,L1,L1,L1;
             L2,L3,L3,L3;
             L3,L2,L3,L3;
             L3,L3,L2,L3;
             L3,L3,L3,L2];
     weight3=[W1;W2;W2;W2;W2];

 elseif ngl==11 % 11-point rule
     % Center Point	
     L1=.250000000000000;
     W1=0.013155555555555;
     % Corner Point	
     L2=.0714285714285714;
     L3=.785714285714286;
     W2=0.007622222222222;
     % Edge Center Points	
     L4=0.399403576166799;
     L5=0.100596423833201;
     W3=0.024888888888888;

     point3=[L1,L1,L1,L1;
             L2,L3,L3,L3;
             L3,L2,L3,L3;
             L3,L3,L2,L3;
             L3,L3,L3,L2;
             L4,L4,L5,L5;
             L5,L4,L4,L5;
             L4,L5,L4,L5;
             L4,L5,L5,L4;
             L5,L4,L5,L4;
             L5,L5,L4,L4];
     weight3=[W1;W2;W2;W2;W2;W3;W3;W3;W3;W3;W3];
 else
     disp("Wrong input for 'ngl' value");
 end
 


