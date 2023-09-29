function [shapet10,dhdrt10,dhdst10,dhdtt10]=feisot10(rvalue,svalue,tvalue)

%------------------------------------------------------------------------
%  Purpose:
%     compute isoparametric ten-node tetrahedra shape functions
%     and their derivatves at the selected (integration) point
%     in terms of the natural coordinate 
%
%  Synopsis:
%     [shapet10,dhdrt10,dhdst10,dhdtt10]=feisot10(rvalue,svalue,tvalue)  
%
%  Variable Description:
%     shapet10 - shape functions for three-node element
%     dhdrt10 - derivatives of the shape functions w.r.t. r
%     dhdst10 - derivatives of the shape functions w.r.t. s
%     dhdtt10 - derivatives of the shape functions w.r.t. t
%     rvalue - r coordinate value of the selected point   
%     svalue - s coordinate value of the selected point
%     tvalue - t coordinate value of the selected point
%------------------------------------------------------------------------
 L1=1-rvalue-svalue-tvalue;
 L2=rvalue;
 L3=svalue;
 L4=tvalue;
% shape functions 
shapet10(1) = L1*(2*L1-1.0);
shapet10(2) = L2*(2*L2-1.0);
shapet10(3) = L3*(2*L3-1.0);
shapet10(4) = L4*(2*L4-1.0);
shapet10(5) = 4.0*L2*L1;
shapet10(6) = 4.0*L2*L3;
shapet10(7) = 4.0*L3*L1;
shapet10(8) = 4.0*L4*L1;
shapet10(9) = 4.0*L2*L4;
shapet10(10)= 4.0*L3*L4;

% derivatives
dhdrt10(1) =-4*L1+1.0;
dhdrt10(2) = 4*L2-1.0;
dhdrt10(3) = 0.0;
dhdrt10(4) = 0.0;
dhdrt10(5) = 4*L1-4*L2;
dhdrt10(6) = 4*L3;
dhdrt10(7) =-4*L3;
dhdrt10(8) =-4*L4;
dhdrt10(9) = 4*L4;
dhdrt10(10)= 0.0;

dhdst10(1) = -4*L1+1.0;
dhdst10(2) = 0.0;
dhdst10(3) = 4*L3-1.0;
dhdst10(4) = 0.0;
dhdst10(5) =-4*L2;
dhdst10(6) = 4*L2;
dhdst10(7) = 4*L1-4*L3;
dhdst10(8) =-4*L4;
dhdst10(9) = 0.0;
dhdst10(10)= 4*L4;

dhdtt10(1) =-4*L1+1.0;
dhdtt10(2) = 0.0;
dhdtt10(3) = 0.0;
dhdtt10(4) = 4*L4-1.0;
dhdtt10(5) = -4*L2;
dhdtt10(6) = 0.0;
dhdtt10(7) =-4*L3;
dhdtt10(8) = 4*L1-4*L4;
dhdtt10(9) = 4*L2;
dhdtt10(10)= 4*L3;
 

