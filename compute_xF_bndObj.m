%get
xF = qqd(end,:)';
qF = xF(1:S1.ndof);
qdF = xF(S1.ndof+1:2*S1.ndof);

g = S1.FwdKinematics(qF);
r_des = g(end-3:end-1,4);
R     = g(end-3:end-1,1:3);
eta = S1.ScrewVelocity(qF,qdF);
v_des = R*eta(end-2:end);


qddu = FDU_IDM(S1,0.5,qF,qdF,uqt);
uF = qddu(end-6:end);