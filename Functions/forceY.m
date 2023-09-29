function fy=forceY(t)
%------------------------------------------------------------------------
%  Purpose:
%     An amplitude to determine the Y component of the force which
%     determined in the question
%
%  Synopsis:
%     fy=forceY(t)
%
%  Variable Description:
%     t - time (pseudo or real)
%------------------------------------------------------------------------
if t<=1
    fy=20*t^3+3*t;
else
    fy=20*(2-t)^3+3*(2-t);
end

return
end
