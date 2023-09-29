function [ff]=feasmbl2(ff,f,index)
%----------------------------------------------------------
%  Purpose:
%     Assembly of element vetors into the system matrix
%
%  Synopsis:
%     [ff]=feasmbl2(ff,f,index)
%
%  Variable Description:
%     ff - system vector
%     f  - element vector
%     index - d.o.f. vector associated with an element
%-----------------------------------------------------------

edof = length(index);
for i=1:edof
    ii=index(i);
    ff(ii)=f(i)+ff(ii);
end