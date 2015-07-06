function dzdt = RicSolit(t,z,d1,d2,A2,A3,e,c)
dzdt = [z(2);
        (-(d1-1)*z(2)^2)/z(1)-(d2*z(2)*z(4))/z(3)+z(2)*z(6)+(d1-1)/z(1)+(A3*z(1)*3)/(d1*z(3)^4)+(e*z(1))/2;
        z(4);
        (-d1*z(2)*z(4))/z(1)-((d2-1)*z(4)^2)/z(3)+z(4)*z(6)+A2/(d2*z(3))-2*(A3*z(1)^2)/(d2*z(3)^3)+(e*z(3))/2;
        z(6);
        -1/3];
        %(-z(6))*((d1*z(2))/(z(1)))+(d2*z(4))/z(3)+z(6)^2+c+e*z(5)];
% Ricci Soliton Equations
%   To carry out our numerical study we shall introduce new variables *z1,
%   z2, z3, z4, z5, z6) = (f1, df1/dt, f2, df2/dt, du, du/dt) to the
%   system.
end