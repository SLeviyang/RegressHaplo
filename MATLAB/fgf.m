% FPG routine for the SL problem 




% Initialization of the primal variables

%y = ones(ndat,1)/ndat;
xp = y;
x = xp;






cnt = 1;



%mu = 0;

Pyvec = 2*Pmat'*yvec;
if rSL==0
    QQf = 2*(Pmat'*Pmat);
else
    QQf = 2*(Pmat'*Pmat+rSL*qpM);
end
yvecyvec = yvec'*yvec;
ydatydat = ydat*ydat';

QQ = QQf + kk*ydatydat;

qq1 = (y>0);
qq2 = y(qq1);
QQs = QQ(:,qq1);
QQy = QQs*qq2;
%QQy = QQ*y;
% EQUAL TO LINE BELOW.  CONVOLUTED COMPUTATION!
grad = -Pyvec + QQy - (mu+kk)*ydat; 
%grad = -2*Pmat'*(yvec-Pmat*y) + 2*rSL*qpM*y - mu*ydat + kk*(ydat'*y-1)*ydat;

grad_save = [grad_save;grad]

chkstop

% rr is c(x) - 1
rr = ydat'*y-1;

% WHERE IS VAL FIRST DEFINED? IN chkstop
recnu = max(val,abs(rr));
%recnu = val
ocnt = 0;
rec(ocnt+1) = recnu;
%rr = 1e4; % |g(alpha)|

%rr1 = 1e4;


dvp = 1e10;


totcnt = 0;

fprintf(1,'Iteration: %i, Accuracy:  %e \n',totcnt, recnu);


% OUTER LOOP, which changes kk and mu, stops when the constraint is within glaveps of
% 0.  And stops when ocnt (outer_count) is 1000.
while ((recnu > glaveps))*(ocnt<1e3)*(ocnt<1000)

t = 1;


QQ = QQf + kk*ydatydat;


% qqp = qq1;
% qq1 = (y>0);
% 
% if qqp==qq1
%     qq2 = y(qq1);
% else
%     QQs = QQ(:,qq1);
%     qq2 = y(qq1);
% end

%same as above
qq1 = (y>0);
qq2 = y(qq1);
QQs = QQ(:,qq1);

QQy = QQs*qq2;
%QQy = QQ*y;
grad = -Pyvec + QQy - (mu+kk)*ydat; 


% Py = Pmat*y;
% qpMy = qpM*y;
 ydaty = ydat'*y;
% grad = -2*Pmat'*(yvec-Py) + 2*rSL*qpMy - mu*ydat + kk*(ydaty-1)*ydat;

chkstop


% mu is lagrange multiplier, kk is penalty term
%dv = (yvec-Py)'*(yvec-Py) + rSL*y'*qpMy - mu*(ydaty-1) + 0.5*kk*((ydaty-1)^2);

dv = 0.5*y'*(-2*Pyvec + QQy) + yvecyvec - mu*(ydaty-1) + kk*(0.5-ydaty);
    
%L = mdat*(kk+0.5); = 
%L = max(abs(eig(2*(Pmat'*Pmat + qpM + kk*ydat*ydat'))));
%L = 3e4;
%L = max(eig(qpdat+kk*ydat'*ydat))

estL
L = 1.25*noq;


cnt = 0;    
ocnt = ocnt + 1;
fprintf('in outerloop %d %f\n', ocnt, y(1))


%INNER LOOP,  cnt = inner count loop, totcnt counts inner and outer loops
%INNER LOOP optimizes x for fixed mu and kk
while ((val > 0.1*recnu)+(totcnt==0))*(cnt < 1e10) % stopping criteria for the inner loop

    
totcnt = totcnt + 1;    
cnt = cnt + 1;
%fprintf('in inner loop %d %f \n', cnt, y(1))

%dv - dopt




x = y - grad/L;
x = max(x,0);
%x = min(x,C');

nt = (1+sqrt(1+4*t^2))/2;

%%%Temporarily disabled

y = x + (x - xp)*(t-1)/nt;

%y = x;

xp = x;

t = nt;



%grad = qpdat*y - c + mu*ydat' + kk*(ydat*y)*ydat';
% Py = Pmat*y;
% qpMy = qpM*y;
 ydaty = ydat'*y;
% grad = -2*Pmat'*(yvec-Py) + 2*rSL*qpMy - mu*ydat + kk*(ydaty-1)*ydat;

%QQ = QQf + kk*ydatydat;

qqp = qq1;
qq1 = (y>0);
if qqp==qq1
    qq2 = y(qq1);
else
    QQs = QQ(:,qq1);
    qq2 = y(qq1);
end

% qq1 = (y>0);
% qq2 = y(qq1);
% QQs = QQ(:,qq1);



QQy = QQs*qq2;

%QQy = QQ*y;


grad = -Pyvec + QQy - (mu+kk)*ydat; 



dvp = dv;
%dv = 0.5*y'*qpdat*y - c'*y + mu*ydat*y + 0.5*kk*(ydat*y)^2;
%dv = (yvec-Py)'*(yvec-Py) + rSL*y'*qpMy - mu*(ydaty-1) + 0.5*kk*((ydaty-1)^2);
dv = 0.5*y'*(-2*Pyvec + QQy) + yvecyvec - mu*(ydaty-1) + kk*(0.5-ydaty);
 

chkstop



end %of inner loop



rr = ydat'*x-1;

% augmented lagrangian multiplier update
mu = mu - kk*rr;
% penalty update
kk = kk*1.01;


totcnt = totcnt + 1;    



recnu = max(val,abs(rr));
%fprintf(1,'Iteration: %i, Accuracy:  %e, Objective value: %e   \n',totcnt, recnu, dv);

rec(ocnt+1) = recnu;
rec1(ocnt) = cnt;
rec2(ocnt+1) = abs(rr);




end

