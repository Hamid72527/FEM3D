function [shapet4,dhdrt4,dhdst4,dhdtt4]=feisot4(rvalue,svalue,tvalue)

%------------------------------------------------------------------------
%  Purpose:
%     compute isoparametric four-node tetrahedra shape functions
%     and their derivatves at the selected (integration) point
%     in terms of the natural coordinate 
%
%  Synopsis:
%     [shapet4,dhdrt4,dhdst4,dhdtt4]=feisot4(rvalue,svalue,tvalue)  
%
%  Variable Description:
%     shapet3 - shape functions for three-node element
%     dhdrt3 - derivatives of the shape functions w.r.t. r
%     dhdst3 - derivatives of the shape functions w.r.t. s
%     dhdtt4 - derivatives of the shape functions w.r.t. t
%     rvalue - r coordinate value of the selected point   
%     svalue - s coordinate value of the selected point
%     tvalue - t coordinate value of the selected point
%  Notes:
%     1st node at (0,0,0), 2nd node at (1,0,0),
%     3rd node at (0,1,0), 4th node at (0,0,1)
%------------------------------------------------------------------------

% shape functions
 shapet4(1)=1-rvalue-svalue-tvalue;
 shapet4(2)=rvalue;
 shapet4(3)=svalue;
 shapet4(4)=tvalue;

% derivatives
 dhdrt4(1)=-1;
 dhdrt4(2)=1;
 dhdrt4(3)=0;
 dhdrt4(4)=0;

 dhdst4(1)=-1;
 dhdst4(2)=0;
 dhdst4(3)=1;
 dhdst4(4)=0;
 
 dhdtt4(1)=-1;
 dhdtt4(2)=0;
 dhdtt4(3)=0;
 dhdtt4(4)=1;
