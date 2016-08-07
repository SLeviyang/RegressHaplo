be = 1e10;
qq = zeros(ndat,1);
for i = 1:ndat
    qqq = qq;
    qqq(i) = 1;
    Pq = Pmat(:,i);
    qpMq = qpM(:,i);
    can = (yvec-Pq)'*(yvec-Pq)+rSL*qpMq(i);
    if (can)<be
        be = can;
        bev = qqq;
        ir = i;
    end
end
        