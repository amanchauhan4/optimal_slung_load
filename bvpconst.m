function res = bvpconst(ya, yb, x_initial, x_final)
global beta
r1 = ya(1:16)-x_initial;
% r2 = (yb(1:16)-x_final);
r2 = yb(17:32) - 2*beta*(yb(1:16) - x_final);
% ustar = UStarFormula(yb(1:16), yb(17:32), constants);
% xdotfinal = bvpode(1, yb, constants);
% r3 =    1 + 0.5*transpose(ustar)*R*ustar + ...
%         0.5*transpose(yb(1:16))*[W1, zeros(8, 8); zeros(8, 8), W2]*yb(1:16) + ...
%         transpose(yb(17:32))*xdotfinal(1:16);
res = [r1;r2];
end