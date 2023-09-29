function [shapet20,dhdrt20,dhdst20,dhdtt20]=feisot20(rvalue,svalue,tvalue)

%------------------------------------------------------------------------
%  Purpose:
%     compute isoparametric twenty-node tetrahedra shape functions
%     and their derivatves at the selected (integration) point
%     in terms of the natural coordinate 
%
%  Synopsis:
%     [shapet20,dhdrt20,dhdst20,dhdtt20]=feisot20(rvalue,svalue,tvalue)  
%
%  Variable Description:
%     shapet20 - shape functions for three-node element
%     dhdrt20 - derivatives of the shape functions w.r.t. r
%     dhdst20 - derivatives of the shape functions w.r.t. s
%     dhdtt20 - derivatives of the shape functions w.r.t. t
%     rvalue - r coordinate value of the selected point   
%     svalue - s coordinate value of the selected point
%     tvalue - t coordinate value of the selected point
%  Notes:
%     1st node at (0,0,0), 2nd node at (1,0,0),
%     3rd node at (0,1,0), 4th node at (0,0,1), ...
%--------------------------------------------------------------------------
 L1=1-rvalue-svalue-tvalue;
 L2=rvalue;
 L3=svalue;
 L4=tvalue;
% shape functions ---------------------------------------------------------
% corner nodes
 shapet20(1)=0.5*(3*L1-1)*(3*L1-2)*L1;
 shapet20(2)=0.5*(3*L2-1)*(3*L2-2)*L2;
 shapet20(3)=0.5*(3*L3-1)*(3*L3-2)*L3;
 shapet20(4)=0.5*(3*L4-1)*(3*L4-2)*L4;
 
% edge nodes
 shapet20(5)=(9/2)*(3*L1-1)*L1*L3;
 shapet20(6)=(9/2)*(3*L3-1)*L1*L3;
 shapet20(7)=(9/2)*(3*L1-1)*L1*L2;
 shapet20(8)=(9/2)*(3*L2-1)*L1*L2;
 shapet20(9)=(9/2)*(3*L2-1)*L2*L3;
 shapet20(10)=(9/2)*(3*L3-1)*L2*L3;
 shapet20(11)=(9/2)*(3*L1-1)*L1*L4;
 shapet20(12)=(9/2)*(3*L4-1)*L1*L4;
 shapet20(13)=(9/2)*(3*L2-1)*L2*L4;
 shapet20(14)=(9/2)*(3*L4-1)*L2*L4;
 shapet20(15)=(9/2)*(3*L3-1)*L3*L4;
 shapet20(16)=(9/2)*(3*L4-1)*L3*L4;

 % centre surface nodes
 shapet20(17)=27*L2*L3*L4;
 shapet20(18)=27*L1*L2*L3;
 shapet20(19)=27*L1*L3*L4;
 shapet20(20)=27*L1*L2*L4;
%--------------------------------------------------------------------------
% derivatives
 dhdrt20(1)=18*L2 + 18*L3 + 18*L4 - (27*L2^2)/2 - (27*L3^2)/2 - (27*L4^2)/2 ...
          - 27*L2*L3 - 27*L2*L4 - 27*L3*L4 - 11/2;
 dhdrt20(2)=(27*L2^2)/2 - 9*L2 + 1;
 dhdrt20(3)=0;
 dhdrt20(4)=0;
 dhdrt20(5)=(9*L3*(6*L2 + 6*L3 + 6*L4 - 5))/2;
 dhdrt20(6)=-L3*((27*L3)/2 - 9/2);
 dhdrt20(7)= (81*L2^2)/2 + 54*L2*L3 + 54*L2*L4 - 45*L2 + (27*L3^2)/2 ... 
            + 27*L3*L4 - (45*L3)/2 + (27*L4^2)/2 - (45*L4)/2 + 9;
 dhdrt20(8)=36*L2 + (9*L3)/2 + (9*L4)/2 - (81*L2^2)/2 - 27*L2*L3 - 27*L2*L4 - 9/2;
 dhdrt20(9)=(9*L3*(6*L2 - 1))/2;
 dhdrt20(10)=L3*((27*L3)/2 - 9/2);
 dhdrt20(11)=(9*L4*(6*L2 + 6*L3 + 6*L4 - 5))/2;
 dhdrt20(12)=-L4*((27*L4)/2 - 9/2);
 dhdrt20(13)=(9*L4*(6*L2 - 1))/2;
 dhdrt20(14)=L4*((27*L4)/2 - 9/2);
 dhdrt20(15)=0;
 dhdrt20(16)=0;
 dhdrt20(17)=27*L3*L4;
 dhdrt20(18)=-27*L3*(2*L2 + L3 + L4 - 1);
 dhdrt20(19)=-27*L3*L4;
 dhdrt20(20)=-27*L4*(2*L2 + L3 + L4 - 1);
% ------------------------------
 dhdst20(1)=18*L2 + 18*L3 + 18*L4 - (27*L2^2)/2 - (27*L3^2)/2 - (27*L4^2)/2 ...
          - 27*L2*L3 - 27*L2*L4 - 27*L3*L4 - 11/2;
 dhdst20(2)=0;
 dhdst20(3)=(27*L3^2)/2 - 9*L3 + 1;
 dhdst20(4)=0;
 dhdst20(5)=(27*L2^2)/2 + 54*L2*L3 + 27*L2*L4 - (45*L2)/2 + (81*L3^2)/2 ...
           + 54*L3*L4 - 45*L3 + (27*L4^2)/2 - (45*L4)/2 + 9;
 dhdst20(6)=(9*L2)/2 + 36*L3 + (9*L4)/2 - (81*L3^2)/2 - 27*L2*L3 - 27*L3*L4 - 9/2;
 dhdst20(7)=(9*L2*(6*L2 + 6*L3 + 6*L4 - 5))/2;
 dhdst20(8)=-L2*((27*L2)/2 - 9/2);
 dhdst20(9)=L2*((27*L2)/2 - 9/2);
 dhdst20(10)=(9*L2*(6*L3 - 1))/2;
 dhdst20(11)=(9*L4*(6*L2 + 6*L3 + 6*L4 - 5))/2;
 dhdst20(12)=-L4*((27*L4)/2 - 9/2);
 dhdst20(13)=0;
 dhdst20(14)=0;
 dhdst20(15)=(9*L4*(6*L3 - 1))/2;
 dhdst20(16)=L4*((27*L4)/2 - 9/2);
 dhdst20(17)=27*L2*L4;
 dhdst20(18)=-27*L2*(L2 + 2*L3 + L4 - 1);
 dhdst20(19)=-27*L4*(L2 + 2*L3 + L4 - 1);
 dhdst20(20)=-27*L2*L4;
% ------------------------------
 dhdtt20(1)=18*L2 + 18*L3 + 18*L4 - (27*L2^2)/2 - (27*L3^2)/2 - (27*L4^2)/2 ...
          - 27*L2*L3 - 27*L2*L4 - 27*L3*L4 - 11/2;
 dhdtt20(2)=0;
 dhdtt20(3)=0;
 dhdtt20(4)=(27*L4^2)/2 - 9*L4 + 1;
 dhdtt20(5)=(9*L3*(6*L2 + 6*L3 + 6*L4 - 5))/2;
 dhdtt20(6)=-L3*((27*L3)/2 - 9/2);
 dhdtt20(7)=(9*L2*(6*L2 + 6*L3 + 6*L4 - 5))/2;
 dhdtt20(8)=-L2*((27*L2)/2 - 9/2);
 dhdtt20(9)=0;
 dhdtt20(10)=0;
 dhdtt20(11)= (27*L2^2)/2 + 27*L2*L3 + 54*L2*L4 - (45*L2)/2 + (27*L3^2)/2 ...
             + 54*L3*L4 - (45*L3)/2 + (81*L4^2)/2 - 45*L4 + 9;
 dhdtt20(12)=(9*L2)/2 + (9*L3)/2 + 36*L4 - (81*L4^2)/2 - 27*L2*L4 - 27*L3*L4 - 9/2;
 dhdtt20(13)=L2*((27*L2)/2 - 9/2);
 dhdtt20(14)=(9*L2*(6*L4 - 1))/2;
 dhdtt20(15)=L3*((27*L3)/2 - 9/2);
 dhdtt20(16)=(9*L3*(6*L4 - 1))/2;
 dhdtt20(17)=27*L2*L3;
 dhdtt20(18)=-27*L2*L3;
 dhdtt20(19)=-27*L3*(L2 + L3 + 2*L4 - 1);
 dhdtt20(20)=-27*L2*(L2 + L3 + 2*L4 - 1);