function fx=forceX(t)
%------------------------------------------------------------------------
%  Purpose:
%     An amplitude to determine the X component of the force which
%     determined in the question
%
%  Synopsis:
%     fx=forceX(t)
%
%  Variable Description:
%     t - time (pseudo or real)
%------------------------------------------------------------------------
if t<=1
    fx=10*t^3+t;
else
    fx=10*(2-t)^3+(2-t);
end

return
end
