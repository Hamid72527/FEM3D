function fz=forceZ(t)
%------------------------------------------------------------------------
%  Purpose:
%     An amplitude to determine the Z component of the force which
%     determined in the question
%
%  Synopsis:
%     fy=forceZ(t)
%
%  Variable Description:
%     t - time (pseudo or real)
%------------------------------------------------------------------------
if t<=1
    fz=(4*t);
else
    fz=(4*(2-t));
end

return
end
