clear

glaveps = 1e-6
fprintf(1,'Reading data from file...   \n'); 

FN1 = strcat('../Igor_P5.csv');
Pmatf=csvread(FN1);

FN2 = strcat('../Igor_y5.csv');

% y VECTOR LOADED
yvec=csvread(FN2);

% P MATRIX DEFINED
Pmat = Pmatf;

[mdat,ndat] = size(Pmat)

ndatst = ndat;



fprintf(1,'Forming matrix M and QP problem...   \n'); 

% M MATRIX CONSTRUCTED
ydat = ones(ndat,1);
qpM = ones(ndat,ndat);
qpM = qpM - diag(ydat);

rSL = 5;
kk = 1;
mu = 1;

FN3 = strcat('../debug_pi.csv');
y=csvread(FN3);


y_save = []
grad_save = []

fgf