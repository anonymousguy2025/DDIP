function [Res,Jac] = IDResidueJacobian(Tr,t,qu,dt,qqdqdd_pre,uqt)

beta = 0.25;
gamma = 0.5;

ndof = Tr.ndof;

q_pre = qqdqdd_pre(1:Tr.ndof);
qd_pre = qqdqdd_pre(Tr.ndof+1:2*Tr.ndof);
qdd_pre = qqdqdd_pre(2*Tr.ndof+1:3*Tr.ndof);

q = qu(1:Tr.ndof);
u = qu(Tr.ndof+1:Tr.ndof+Tr.nact);

qd  = gamma/(beta*dt)*(q-q_pre)+(1-gamma/beta)*qd_pre+dt*(1-gamma/(2*beta))*qdd_pre;
qdd = 1/(beta*dt^2)*(q-q_pre)-1/(beta*dt)*qd_pre+(1-1/(2*beta))*qdd_pre;

N    = Tr.N;
nsig = Tr.nsig;
nj   = nsig-N; %since there are 3 soft links and 1 rigid link: 3*(nip-1)+3 joints for soft link, 1 joint for rigid link

iLpre = Tr.iLpre; %previous link
g_ini = Tr.g_ini; %g_ini of each link
%% Forward pass, evaluate joint quantities and kinematic, differential kinematic and derivative terms

%Joint quantities
h = zeros(1,nj);
Omega = zeros(6,nj);
Z = zeros(6*nj,ndof);
gstep = zeros(4*nj,4);
Adgstepinv = zeros(6*nj,6);
T = zeros(6*nj,6);
S = zeros(6*nj,ndof);
Sd = zeros(6*nj,ndof);
f = zeros(4,nj);
fd = zeros(4,nj);
adjOmegap = zeros(24*nj,6); %Powers of adjOmega, used later
dSdq_qd = zeros(6*nj,ndof);
dSdq_qdd = zeros(6*nj,ndof);
dSddq_qd = zeros(6*nj,ndof);
R = zeros(6*nj,ndof);
Q = zeros(6*nj,ndof);
Y = zeros(6*nj,ndof);

g   = zeros(4*nsig,4);
eta = zeros(6*nsig,1); %J*qd
etadot = zeros(6*nsig,1); %Jd*qd

S_B = zeros(6*nsig,ndof); %Jacobian
R_B = zeros(6*nsig,ndof);
Q_B = zeros(6*nsig,ndof);
Y_B = zeros(6*nsig,ndof);

g_tip  = repmat(eye(4),N,1);
S_Btip = repmat(zeros(6,ndof),N,1);
R_Btip = repmat(zeros(6,ndof),N,1);
Q_Btip = repmat(zeros(6,ndof),N,1);
Y_Btip = repmat(zeros(6,ndof),N,1);
eta_tip = zeros(N*6,1);
etadot_tip = zeros(N*6,1);

j = 1; % only 1 soft division
ij = 1;
isig = 1;
dof_start = 1;

for i=1:N

    if iLpre(i)>0
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        S_Bhere      = Ad_g_ini_inv*S_Btip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        R_Bhere      = Ad_g_ini_inv*R_Btip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        Q_Bhere      = Ad_g_ini_inv*Q_Btip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        Y_Bhere      = Ad_g_ini_inv*Y_Btip((iLpre(i)-1)*6+1:iLpre(i)*6,:);

        g_here       = g_tip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        eta_here     = Ad_g_ini_inv*eta_tip((iLpre(i)-1)*6+1:iLpre(i)*6);
        etadot_here  = Ad_g_ini_inv*etadot_tip((iLpre(i)-1)*6+1:iLpre(i)*6);
    else 
        S_Bhere  = zeros(6,ndof);
        R_Bhere  = zeros(6,ndof);
        Q_Bhere  = zeros(6,ndof);
        Y_Bhere  = zeros(6,ndof);

        g_here = g_ini((i-1)*4+1:i*4,:);
        eta_here = zeros(6,1);
        etadot_here = zeros(6,1);
    end

    %0 of joint
    
    S_B((isig-1)*6+1:6*isig,:) = S_Bhere;
    R_B((isig-1)*6+1:6*isig,:) = R_Bhere;
    Q_B((isig-1)*6+1:6*isig,:) = Q_Bhere;
    Y_B((isig-1)*6+1:6*isig,:) = Y_Bhere;

    g((isig-1)*4+1:4*isig,:) = g_here;
    eta((isig-1)*6+1:6*isig) = eta_here;
    etadot((isig-1)*6+1:6*isig) = etadot_here;

    dof_here = Tr.CVTwists{i}(1).dof;
    dofs_here = dof_start:dof_start+dof_here-1;
    Phi = Tr.CVTwists{i}(1).B;
    h(ij) = 1;
    
    [Omega(:,ij),Z((ij-1)*6+1:6*ij,dofs_here),gstep((ij-1)*4+1:4*ij,:),T((ij-1)*6+1:6*ij,:),S((ij-1)*6+1:6*ij,dofs_here),Sd((ij-1)*6+1:6*ij,dofs_here),...
        f(:,ij),fd(:,ij),adjOmegap((ij-1)*24+1:24*ij,:),dSdq_qd((ij-1)*6+1:6*ij,dofs_here),dSdq_qdd((ij-1)*6+1:6*ij,dofs_here),dSddq_qd((ij-1)*6+1:6*ij,dofs_here)]...
        =JointQuantitiesR_mex(Phi,dof_here,q(dofs_here),qd(dofs_here),qdd(dofs_here));
    Adgstepinv((ij-1)*6+1:6*ij,:) = dinamico_Adjoint(ginv(gstep((ij-1)*4+1:4*ij,:)));

    eta_plus_here    = eta_here+S((ij-1)*6+1:6*ij,dofs_here)*qd(dofs_here);
    etadot_plus_here = etadot_here+S((ij-1)*6+1:6*ij,dofs_here)*qdd(dofs_here)+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,dofs_here)*qd(dofs_here)+Sd((ij-1)*6+1:6*ij,dofs_here)*qd(dofs_here);

    R((ij-1)*6+1:6*ij,dofs_here) = dinamico_adj(eta_plus_here)*S((ij-1)*6+1:6*ij,dofs_here)+dSdq_qd((ij-1)*6+1:6*ij,dofs_here);
    Q((ij-1)*6+1:6*ij,dofs_here) = dinamico_adj(etadot_plus_here-dinamico_Adjoint(ginv(g_here))*Tr.G)*S((ij-1)*6+1:6*ij,dofs_here)...
           +dinamico_adj(eta_plus_here)*R((ij-1)*6+1:6*ij,dofs_here)+dinamico_adj(eta_here)*dSdq_qd((ij-1)*6+1:6*ij,dofs_here)...
           +dSddq_qd((ij-1)*6+1:6*ij,dofs_here)+dSdq_qdd((ij-1)*6+1:6*ij,dofs_here);
    Y((ij-1)*6+1:6*ij,dofs_here) = R((ij-1)*6+1:6*ij,dofs_here)+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,dofs_here)+Sd((ij-1)*6+1:6*ij,dofs_here);

    % from 0 to 1 of rigid joint
    S_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(S_Bhere+S((ij-1)*6+1:6*ij,:));
    R_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(R_Bhere+R((ij-1)*6+1:6*ij,:));
    Q_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(Q_Bhere+Q((ij-1)*6+1:6*ij,:));
    Y_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(Y_Bhere+Y((ij-1)*6+1:6*ij,:));

    g_here = g_here*gstep((ij-1)*4+1:4*ij,:);
    eta_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(eta_plus_here);
    etadot_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(etadot_plus_here);

    ij = ij+1;
    isig = isig+1;
    dof_start = dof_start+dof_here;

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        %joint to CM
        gi        = Tr.VLinks(Tr.LinkIndex(i)).gi; 
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));

        S_Bhere = Ad_gi_inv*S_Bhere;
        R_Bhere = Ad_gi_inv*R_Bhere;
        Q_Bhere = Ad_gi_inv*Q_Bhere;
        Y_Bhere = Ad_gi_inv*Y_Bhere;

        g_here = g_here*gi;
        eta_here = Ad_gi_inv*eta_here;
        etadot_here = Ad_gi_inv*etadot_here;

        S_B((isig-1)*6+1:6*isig,:) = S_Bhere;
        R_B((isig-1)*6+1:6*isig,:) = R_Bhere;
        Q_B((isig-1)*6+1:6*isig,:) = Q_Bhere;
        Y_B((isig-1)*6+1:6*isig,:) = Y_Bhere;
        
        g((isig-1)*4+1:4*isig,:) = g_here;
        eta((isig-1)*6+1:6*isig) = eta_here;
        etadot((isig-1)*6+1:6*isig) = etadot_here;

        isig = isig+1;

        %CM to tip
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf; 
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));

        S_Bhere = Ad_gf_inv*S_Bhere;
        R_Bhere = Ad_gf_inv*R_Bhere;
        Q_Bhere = Ad_gf_inv*Q_Bhere;
        Y_Bhere = Ad_gf_inv*Y_Bhere;

        g_here = g_here*gf;
        eta_here = Ad_gf_inv*eta_here;
        etadot_here = Ad_gf_inv*etadot_here;
        
    else
        %gi and gf of soft link are identity hence no transformation is done here
        %at X=0
        S_B((isig-1)*6+1:6*isig,:) = S_Bhere;
        R_B((isig-1)*6+1:6*isig,:) = R_Bhere;
        Q_B((isig-1)*6+1:6*isig,:) = Q_Bhere;
        Y_B((isig-1)*6+1:6*isig,:) = Y_Bhere;
        
        g((isig-1)*4+1:4*isig,:) = g_here;
        eta((isig-1)*6+1:6*isig) = eta_here; % to find FD
        etadot((isig-1)*6+1:6*isig) = etadot_here;
        
        isig = isig+1;
        
        dof_here = Tr.CVTwists{i}(j+1).dof;
        dofs_here = dof_start:dof_start+dof_here-1;
        nip = Tr.CVTwists{i}(j+1).nip;
        Lscale = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        h(ij:ij+nip-2) = (Tr.CVTwists{i}(j+1).Xs(2:end)-Tr.CVTwists{i}(j+1).Xs(1:end-1))*Lscale;

        for ii=1:nip-1
            Phi_Z1 = Tr.CVTwists{i}(j+1).B_Z1(6*(ii-1)+1:6*ii,:);
            Phi_Z2 = Tr.CVTwists{i}(j+1).B_Z2(6*(ii-1)+1:6*ii,:);
            xi_star_Z1 = Tr.CVTwists{i}(j+1).xi_star(6*(ii-1)+1:6*ii,2); 
            xi_star_Z2 = Tr.CVTwists{i}(j+1).xi_star(6*(ii-1)+1:6*ii,3); 

            [Omega(:,ij),Z((ij-1)*6+1:6*ij,dofs_here),gstep((ij-1)*4+1:4*ij,:),T((ij-1)*6+1:6*ij,:),S((ij-1)*6+1:6*ij,dofs_here),Sd((ij-1)*6+1:6*ij,dofs_here),...
                f(:,ij),fd(:,ij),adjOmegap((ij-1)*24+1:24*ij,:),dSdq_qd((ij-1)*6+1:6*ij,dofs_here),dSdq_qdd((ij-1)*6+1:6*ij,dofs_here),dSddq_qd((ij-1)*6+1:6*ij,dofs_here)]...
                =JointQuantities_mex(h(ij),Phi_Z1,Phi_Z2,xi_star_Z1,xi_star_Z2,q(dofs_here),qd(dofs_here),qdd(dofs_here));
            Adgstepinv((ij-1)*6+1:6*ij,:)  = dinamico_Adjoint(ginv(gstep((ij-1)*4+1:4*ij,:)));

            eta_plus_here    = eta_here+S((ij-1)*6+1:6*ij,dofs_here)*qd(dofs_here);
            etadot_plus_here = etadot_here+S((ij-1)*6+1:6*ij,dofs_here)*qdd(dofs_here)+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,dofs_here)*qd(dofs_here)+Sd((ij-1)*6+1:6*ij,dofs_here)*qd(dofs_here);

            R((ij-1)*6+1:6*ij,dofs_here) = dinamico_adj(eta_plus_here)*S((ij-1)*6+1:6*ij,dofs_here)+dSdq_qd((ij-1)*6+1:6*ij,dofs_here);
            Q((ij-1)*6+1:6*ij,dofs_here) = dinamico_adj(etadot_plus_here-dinamico_Adjoint(ginv(g_here))*Tr.G)*S((ij-1)*6+1:6*ij,dofs_here)...
                   +dinamico_adj(eta_plus_here)*R((ij-1)*6+1:6*ij,dofs_here)+dinamico_adj(eta_here)*dSdq_qd((ij-1)*6+1:6*ij,dofs_here)...
                   +dSddq_qd((ij-1)*6+1:6*ij,dofs_here)+dSdq_qdd((ij-1)*6+1:6*ij,dofs_here);
            Y((ij-1)*6+1:6*ij,dofs_here) = R((ij-1)*6+1:6*ij,dofs_here)+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,dofs_here)+Sd((ij-1)*6+1:6*ij,dofs_here);
            
            S_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(S_Bhere+S((ij-1)*6+1:6*ij,:));
            R_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(R_Bhere+R((ij-1)*6+1:6*ij,:));
            Q_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(Q_Bhere+Q((ij-1)*6+1:6*ij,:));
            Y_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(Y_Bhere+Y((ij-1)*6+1:6*ij,:));
            
            g_here = g_here*gstep((ij-1)*4+1:4*ij,:);
            eta_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(eta_plus_here);
            etadot_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(etadot_plus_here);
            
            S_B((isig-1)*6+1:6*isig,:) = S_Bhere;
            R_B((isig-1)*6+1:6*isig,:) = R_Bhere;
            Q_B((isig-1)*6+1:6*isig,:) = Q_Bhere;
            Y_B((isig-1)*6+1:6*isig,:) = Y_Bhere;
            
            g((isig-1)*4+1:4*isig,:) = g_here;
            eta((isig-1)*6+1:6*isig) = eta_here; % to find FD
            etadot((isig-1)*6+1:6*isig) = etadot_here;

            ij = ij+1;
            isig = isig+1;
        end
        dof_start = dof_start+dof_here;
    end
    % saving tip values
    S_Btip((i-1)*6+1:i*6,:) = S_Bhere;
    R_Btip((i-1)*6+1:i*6,:) = R_Bhere;
    Q_Btip((i-1)*6+1:i*6,:) = Q_Bhere;
    Y_Btip((i-1)*6+1:i*6,:) = Y_Bhere;

    g_tip((i-1)*4+1:i*4,:) = g_here;
    eta_tip((i-1)*6+1:i*6) = eta_here;
    etadot_tip((i-1)*6+1:i*6) = etadot_here;
end

%% Backward Pass, to compute derivatives

isig = nsig;
ij = nj;
dof_start = ndof;

ID = zeros(ndof,1);
dID_dq = zeros(ndof,ndof);
dID_dqd = zeros(ndof,ndof);
dID_dqdd = zeros(ndof,ndof); %M

M_C_tip = repmat(zeros(6,6),N,1);
N_C_tip = repmat(zeros(6,6),N,1);
F_C_tip = repmat(zeros(6,1),N,1);
P_S_tip = repmat(zeros(6,ndof),N,1);
U_S_tip = repmat(zeros(6,ndof),N,1);
V_S_tip = repmat(zeros(6,ndof),N,1);
W_S_tip = repmat(zeros(6,ndof),N,1);

for i=N:-1:1 %backwards

    M_C = M_C_tip(6*(i-1)+1:6*i,:);
    N_C = N_C_tip(6*(i-1)+1:6*i,:);
    F_C = F_C_tip(6*(i-1)+1:6*i);
    P_S = P_S_tip(6*(i-1)+1:6*i,:);
    U_S = U_S_tip(6*(i-1)+1:6*i,:);
    V_S = V_S_tip(6*(i-1)+1:6*i,:);
    W_S = W_S_tip(6*(i-1)+1:6*i,:);

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        %tip to CM
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf; 
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        coAd_gf   = Ad_gf_inv';
        
        M_C = coAd_gf*M_C*Ad_gf_inv;
        N_C = coAd_gf*N_C*Ad_gf_inv;
        F_C = coAd_gf*F_C;
        P_S = coAd_gf*P_S;
        U_S = coAd_gf*U_S;
        V_S = coAd_gf*V_S;
        W_S = coAd_gf*W_S;
        isig = isig-1;

        %CM to base
        gi        = Tr.VLinks(Tr.LinkIndex(i)).gi; 
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        coAd_gi   = Ad_gi_inv';
        
        M_ap1 = Tr.VLinks(Tr.LinkIndex(i)).M; %at alpha + 1.
        N_ap1 = dinamico_coadjbar(M_ap1*eta(isig*6+1:6*(isig+1)))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1-M_ap1*dinamico_adj(eta(isig*6+1:6*(isig+1)));
        F_ap1 = -M_ap1*dinamico_Adjoint(ginv(g(isig*4+1:4*(isig+1),:)))*Tr.G+M_ap1*etadot(isig*6+1:6*(isig+1))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1*eta(isig*6+1:6*(isig+1));

        M_C = coAd_gi*(M_ap1+M_C)*Ad_gi_inv;
        N_C = coAd_gi*(N_ap1+N_C)*Ad_gi_inv;
        F_C = coAd_gi*(F_ap1+F_C);
        P_S = coAd_gi*P_S;
        U_S = coAd_gi*U_S;
        V_S = coAd_gi*V_S;
        W_S = coAd_gi*W_S;

        %for rigid joint
        M_ap1 = zeros(6,6); %at base
        N_ap1 = zeros(6,6);
        F_ap1 = zeros(6,1);
    else
        %assuming gf and gi of soft are identities
        dof_here = Tr.CVTwists{i}(j+1).dof;
        dofs_here = dof_start-dof_here+1:dof_start;
        Ws = Tr.CVTwists{i}(j+1).Ws; %L*w_gauss
        nip = Tr.CVTwists{i}(j+1).nip;
        isig = isig-1;
        for ii=nip-1:-1:1
        
            coAdgstep = Adgstepinv((ij-1)*6+1:6*ij,:)';
        
            M_ap1 = Ws(ii+1)*Tr.CVTwists{i}(2).Ms(ii*6+1:6*(ii+1),:); %at alpha + 1. multiplied with weight
            N_ap1 = dinamico_coadjbar(M_ap1*eta(isig*6+1:6*(isig+1)))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1-M_ap1*dinamico_adj(eta(isig*6+1:6*(isig+1)));
            F_ap1 = -M_ap1*dinamico_Adjoint(ginv(g(isig*4+1:4*(isig+1),:)))*Tr.G+M_ap1*etadot(isig*6+1:6*(isig+1))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1*eta(isig*6+1:6*(isig+1));

            M_C = coAdgstep*(M_ap1+M_C)*Adgstepinv((ij-1)*6+1:6*ij,:);
            N_C = coAdgstep*(N_ap1+N_C)*Adgstepinv((ij-1)*6+1:6*ij,:);
            F_C = coAdgstep*(F_ap1+F_C);

            %split to two steps to optimize computation
            P_S = coAdgstep*P_S;
            U_S = coAdgstep*U_S;
            V_S = coAdgstep*V_S;
            W_S = coAdgstep*W_S;

            P_S(:,dofs_here) = P_S(:,dofs_here)+dinamico_coadjbar(F_C)*S((ij-1)*6+1:6*ij,dofs_here);
            U_S(:,dofs_here) = U_S(:,dofs_here)+N_C*R((ij-1)*6+1:6*ij,dofs_here)+M_C*Q((ij-1)*6+1:6*ij,dofs_here);
            V_S(:,dofs_here) = V_S(:,dofs_here)+N_C*S((ij-1)*6+1:6*ij,dofs_here)+M_C*Y((ij-1)*6+1:6*ij,dofs_here);
            W_S(:,dofs_here) = W_S(:,dofs_here)+M_C*S((ij-1)*6+1:6*ij,dofs_here);

            Phi_Z1 = Tr.CVTwists{i}(j+1).B_Z1(6*(ii-1)+1:6*ii,:);
            Phi_Z2 = Tr.CVTwists{i}(j+1).B_Z2(6*(ii-1)+1:6*ii,:);

            dSTdq_FC = compute_dSTdq_FC_mex(h(ij),Omega(:,ij),Phi_Z1,Phi_Z2,Z((ij-1)*6+1:6*ij,dofs_here),T((ij-1)*6+1:6*ij,:),f(:,ij),fd(:,ij),adjOmegap((ij-1)*24+1:24*ij,:),F_C);
            dSTdq_FC_full = zeros(dof_here,ndof);
            dSTdq_FC_full(:,dofs_here) = dSTdq_FC;
            
            ID(dofs_here,:) = ID(dofs_here,:)+S((ij-1)*6+1:6*ij,dofs_here)'*F_C;
            dID_dq(dofs_here,:)  = dID_dq(dofs_here,:)+dSTdq_FC_full+S((ij-1)*6+1:6*ij,dofs_here)'*(N_C*R_B((isig-1)*6+1:6*isig,:)+M_C*Q_B((isig-1)*6+1:6*isig,:)+U_S+P_S);
            dID_dqd(dofs_here,:) = dID_dqd(dofs_here,:)+S((ij-1)*6+1:6*ij,dofs_here)'*(N_C*S_B((isig-1)*6+1:6*isig,:)+M_C*Y_B((isig-1)*6+1:6*isig,:)+V_S);
            dID_dqdd(dofs_here,:)       = dID_dqdd(dofs_here,:)+S((ij-1)*6+1:6*ij,dofs_here)'*(M_C*S_B((isig-1)*6+1:6*isig,:)+W_S);

            isig = isig-1;
            ij = ij-1;
        
        end
        %for rigid joint
        M_ap1 = Ws(1)*Tr.CVTwists{i}(2).Ms(1:6,:); %at X = 0 usually Ws is 0
        N_ap1 = dinamico_coadjbar(M_ap1*eta(isig*6+1:6*(isig+1)))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1-M_ap1*dinamico_adj(eta(isig*6+1:6*(isig+1)));
        F_ap1 = -M_ap1*dinamico_Adjoint(ginv(g(isig*4+1:4*(isig+1),:)))*Tr.G+M_ap1*etadot(isig*6+1:6*(isig+1))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1*eta(isig*6+1:6*(isig+1));
        dof_start=dof_start-dof_here;
    end
    %joint
    dof_here = Tr.CVTwists{i}(1).dof;

    if dof_here>0

        dofs_here = dof_start-dof_here+1:dof_start;
        coAdgstep = Adgstepinv((ij-1)*6+1:6*ij,:)';
        
        M_C = coAdgstep*(M_ap1+M_C)*Adgstepinv((ij-1)*6+1:6*ij,:);
        N_C = coAdgstep*(N_ap1+N_C)*Adgstepinv((ij-1)*6+1:6*ij,:);
        F_C = coAdgstep*(F_ap1+F_C);
    
        %split to two steps to optimize computation
        P_S = coAdgstep*P_S;
        U_S = coAdgstep*U_S;
        V_S = coAdgstep*V_S;
        W_S = coAdgstep*W_S;
    
        P_S(:,dofs_here) = P_S(:,dofs_here)+dinamico_coadjbar(F_C)*S((ij-1)*6+1:6*ij,dofs_here);
        U_S(:,dofs_here) = U_S(:,dofs_here)+N_C*R((ij-1)*6+1:6*ij,dofs_here)+M_C*Q((ij-1)*6+1:6*ij,dofs_here);
        V_S(:,dofs_here) = V_S(:,dofs_here)+N_C*S((ij-1)*6+1:6*ij,dofs_here)+M_C*Y((ij-1)*6+1:6*ij,dofs_here);
        W_S(:,dofs_here) = W_S(:,dofs_here)+M_C*S((ij-1)*6+1:6*ij,dofs_here);
    
        dSTdq_FC = compute_dSTdq_FCR_mex(dof_here,Omega(:,ij),Z((ij-1)*6+1:6*ij,dofs_here),f(:,ij),fd(:,ij),adjOmegap((ij-1)*24+1:24*ij,:),F_C);
        dSTdq_FC_full = zeros(dof_here,ndof);
        dSTdq_FC_full(:,dofs_here) = dSTdq_FC;
        
        ID(dofs_here,:) = ID(dofs_here,:)+S((ij-1)*6+1:6*ij,dofs_here)'*F_C;
        dID_dq(dofs_here,:)  = dID_dq(dofs_here,:)+dSTdq_FC_full+S((ij-1)*6+1:6*ij,dofs_here)'*(N_C*R_B((isig-1)*6+1:6*isig,:)+M_C*Q_B((isig-1)*6+1:6*isig,:)+U_S+P_S);
        dID_dqd(dofs_here,:) = dID_dqd(dofs_here,:)+S((ij-1)*6+1:6*ij,dofs_here)'*(N_C*S_B((isig-1)*6+1:6*isig,:)+M_C*Y_B((isig-1)*6+1:6*isig,:)+V_S);
        dID_dqdd(dofs_here,:)       = dID_dqdd(dofs_here,:)+S((ij-1)*6+1:6*ij,dofs_here)'*(M_C*S_B((isig-1)*6+1:6*isig,:)+W_S);

    end
    
    isig = isig-1;
    ij = ij-1;
    dof_start=dof_start-dof_here;

    % projecting to the tip of previous link
    ip = iLpre(i);

    if ip>0
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        coAd_g_ini = Ad_g_ini_inv';
    
        M_C_tip(6*(ip-1)+1:6*ip,:) = M_C_tip(6*(ip-1)+1:6*ip,:)+coAd_g_ini*M_C*Ad_g_ini_inv;
        N_C_tip(6*(ip-1)+1:6*ip,:) = N_C_tip(6*(ip-1)+1:6*ip,:)+coAd_g_ini*N_C*Ad_g_ini_inv;
        F_C_tip(6*(ip-1)+1:6*ip)   = F_C_tip(6*(ip-1)+1:6*ip)+coAd_g_ini*F_C;
        P_S_tip(6*(ip-1)+1:6*ip,:) = P_S_tip(6*(ip-1)+1:6*ip,:)+coAd_g_ini*P_S;
        U_S_tip(6*(ip-1)+1:6*ip,:) = U_S_tip(6*(ip-1)+1:6*ip,:)+coAd_g_ini*U_S;
        V_S_tip(6*(ip-1)+1:6*ip,:) = V_S_tip(6*(ip-1)+1:6*ip,:)+coAd_g_ini*V_S;
        W_S_tip(6*(ip-1)+1:6*ip,:) = W_S_tip(6*(ip-1)+1:6*ip,:)+coAd_g_ini*W_S;
    end

end

%% derivatives of internal wrenches
dtau_dq  = -Tr.K;
dtau_dqd = -Tr.D;

B = Tr.Bqj1; %only 1 DOF joints, its derivative is 0

tau = B*u-Tr.K*q-Tr.D*qd;
%% Putting together
qqdqddjt = zeros(Tr.nact,1);
qqdqddjt(1:2) = uqt;

Res = [tau-ID;B'*q-qqdqddjt(:,1)];
Jac = [dtau_dq+dtau_dqd*gamma/(beta*dt)-(dID_dq+dID_dqd*gamma/(beta*dt)+dID_dqdd*1/(beta*dt^2)) B;B' zeros(Tr.nact,Tr.nact)];
end