clear
% main controling module

glaveps = 1e-4;



% Reading the file from a general data format


fprintf(1,'Reading data from file...   \n'); 

FN1 = strcat('Igor_P1.csv');
Pmatf=xlsread(FN1);

FN2 = strcat('Igor_y1.csv');

# y VECTOR LOADED
yvec=xlsread(FN2);

% P MATRIX DEFINED
Pmat = Pmatf;

[mdat,ndat] = size(Pmat)

ndatst = ndat;



fprintf(1,'Forming matrix M and QP problem...   \n'); 

% M MATRIX CONSTRUCTED
ydat = ones(ndat,1);
qpM = ones(ndat,ndat);
qpM = qpM - diag(ydat);


% Selection of k

% k0 = 16;
% kk = k0;

% fprintf(1,'Estimating L...'); 
% egg = eig(2*(Pmat'*Pmat + qpM) + kk*ydat*ydat');
% L = 1.5*max(abs(egg))
% fprintf(1,' Estimated\n'); 

kar = [];

% parameter rho

rhoarray = [0.1];

qstat = [];

sch = 1;

rSL = rhoarray;
    
for k0 = [2 4 8 16]

for trial = 1:100

    
kk = k0;    
    
time1 = cputime;
%%%%%%%%% Presolve if needed

%y = ones(ndat,1)/ndat;
if trial==1
    ndat = ndatst;
    fprintf(1,'Presolving the QP...\n'); 
    quickfind
    y = bev;
    dv = be;
    mu = 0;
elseif trial==2
    ndat = ndatst;
    mu = 0;
    y = ones(ndat,1)/ndat;
elseif (trial>=3) && (trial<=10)
    ndat = ndatst;
    mu = 0;
    y = rand(ndat,1);
else
    ndat = ndatst;
    mu = 0;
    y = rand(ndat,1);
    y = y/norm(y);
end


%%%%%%%%%%% Main subroutine
fprintf(1,'Solving the QP...\n'); 

fgf

fprintf(1,' Solved\n'); 

time2 = cputime;

fprintf(1,'Run number: %i \n',sch);
fprintf(1,'Infeasibility:  %e \n',abs(rr));
fprintf(1,'Accuracy:  %e \n',val);
fprintf(1,'Iterations:  %i \n',totcnt);
fprintf(1,'Solution time:  %e sec\n',time2-time1);
fprintf(1,'\n')

qstat(:,sch) =qq1;
sch = sch + 1;
kar = [kar k0];
end
end

rhost = rSL*ones(1,(sch-1));

%%%%% Running accross different indices
% FOR EVERY SOLUTION pi, COMPUTES THE p SOLUTION WITH NO
% PENALTY
store=x;
storev = dv;

kk = 8;
glaveps = 1e-6;

[qh,qnum] = size(qstat);
# setting penalty to zero
rSL = 0;

xrec = [];
frec = [];

qq1p = (qstat(:,1)>0);

Pmat = Pmatf(:,qq1p);
[mdat,ndat] = size(Pmat)

ydat = ones(ndat,1);
y = ydat/ndat;

% fprintf(1,'Presolving the QP...\n'); 
% quickfind
% y = bev;
% dv = be;
% %mu = grad(ir);

mu = 0;


fgf;
hii = 1;
ndat = length(qq1p);
xs = zeros(ndat,1);
for hi = 1:ndat
    if qq1p(hi)>0
        xs(hi) = x(hii);
        hii = hii + 1;
    end
end
        
xrec = [xrec xs];
frec = [frec dv];
for qi = 2:qnum
    %if max(qstat(:,qi)~=qstat(:,qi-1))
        qq1p = (qstat(:,qi)>0);

        Pmat = Pmatf(:,qq1p);
        [mdat,ndat] = size(Pmat)

        ydat = ones(ndat,1);
        y = ydat/ndat;
%         fprintf(1,'Presolving the QP...\n'); 
%         quickfind
%         y = bev;
%         dv = be;
%mu = grad(ir);
        mu = 0;
        kk = 8;

        
        fgf;
        hii = 1;
        ndat = length(qq1p);
        xs = zeros(ndat,1);
        for hi = 1:ndat
            if qq1p(hi)>0
                xs(hi) = x(hii);
                hii = hii + 1;
            end
        end

        
        xrec = [xrec xs];
        frec = [frec dv];
    %end
end

wxrec = [frec;rhost;sum(qstat);kar;xrec];
csvwrite('ps.csv',wxrec);
