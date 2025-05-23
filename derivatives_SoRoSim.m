function qdqdd = derivatives_SoRoSim(Tr,t,qqd,uqt) %u is a function of t
    qd  = qqd(Tr.ndof+1:2*Tr.ndof);
    qdd = FD_SoRoSim(Tr,t,qqd,uqt);
    qqdqddjt = uqt(t);
    qd(Tr.i_jactq) = qqdqddjt(:,2);
    qdqdd = [qd;qdd];
    disp(t);
end