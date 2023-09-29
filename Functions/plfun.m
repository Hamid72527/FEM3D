function [eSigmay,hmodule]=plfun(ePEEQ)
%------------------------------------------------------------------------
%  Purpose:
%     to determine the yield stress and plastic tangent modulus (H) in terms
%     of plastic accumulated strain
%
%  Synopsis:
%     [eSigmay,hmodule]=plfun(ePEEQ)
%
%  Variable Description:
%     eSigmay - yield stress
%     hmodule - plastic tangent modulus
%     ePEEQ - plastic accumulated strain
%------------------------------------------------------------------------
a =   4.235e+07;
b =       12.15;
c =  -2.235e+07;
d =       -2054;

eSigmay = a*exp(b*ePEEQ) + c*exp(d*ePEEQ);
hmodule=a*b*exp(b*ePEEQ) + c*d*exp(d*ePEEQ);
return
end

