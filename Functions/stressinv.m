function [I1,J2,J3]=stressinv(stress)
%------------------------------------------------------------------------
%  Purpose:
%     Calculate Stress Invarients for two-dimensional and two-dimensional
%     Stress tensord
%
%  Synopsis:
%     [I1,J2,J3]=StressInvCal(stress) 
%
%  Variable Description:
%     stress_dev - the deviatoric stress tensor
%------------------------------------------------------------------------
if length (stress) == 6
    stress=[stress(1) stress(4) stress(6);
            stress(4) stress(2) stress(5);
            stress(6) stress(5) stress(3)];
elseif length (stress) == 9
    stress=[stress(1) stress(4) stress(9);
            stress(5) stress(2) stress(6);
            stress(8) stress(7) stress(3)];
end
    
I1=trace(stress);
stress_dev=stress-(I1/3)*diag(ones(3,1));
J2=0;
for i=1:3
    for j=1:3
        J2=J2+0.5*stress_dev(i,j)*stress_dev(i,j);
    end
end
J3=det(stress_dev);