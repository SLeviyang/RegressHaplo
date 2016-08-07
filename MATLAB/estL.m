%QQ = 2*(Pmat'*Pmat + rSL*qpM) + kk*ydat*ydat';

% q is set to normalized 1 vector
q = ones(ndat,1);
noq = norm(q);
q = q/noq;
i = 0;
noqp = 1e10;
while abs(noq-noqp)/noq > 0.1
    q = QQ*q;
    noqp = noq;
    noq = norm(q);
    q = q/noq;
    i = i + 1;
end
