function qddu = FDU_IDM(Tr,t,q,qd,uqt)
N    = Tr.N;
ndof = Tr.ndof;
nsig = Tr.nsig; %nip+1 for soft link (rigid joint+soft body) and 2 for rigid link
nj   = nsig-N; %nip joints for soft link (rigid joint+distributed joint) and 1 for rigid link


qqdqddjt = uqt(t);
for i=1:Tr.n_jact
    if ~Tr.WrenchControlled(i)
        q(Tr.i_jactq(i))  = qqdqddjt(i,1); %replacing q_k and qdot_k
        qd(Tr.i_jactq(i)) = qqdqddjt(i,2);
    end
end

%% Evaluate Quantities at every joints
gstep = zeros(4*nj,4);
Adgstepinv = zeros(6*nj,6);
S = zeros(6*nj,ndof);
Sd = zeros(6*nj,ndof);

j = 1; % only 1 soft division
ij = 1;
isig = 1;

for i=1:N
    Phi = Tr.PhiE(6*(isig-1)+1:6*isig,:);
    [gstep((ij-1)*4+1:4*ij,:),S((ij-1)*6+1:6*ij,:),Sd((ij-1)*6+1:6*ij,:)]=JointQuantitiesR_FD_mex(Phi,q,qd);
    Adgstepinv((ij-1)*6+1:6*ij,:)  = dinamico_Adjoint(ginv(gstep((ij-1)*4+1:4*ij,:)));
    ij = ij+1;
    isig = isig+1;
    %for soft
    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='s'

        nip = Tr.CVTwists{i}(j+1).nip;
        Lscale = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        h = (Tr.CVTwists{i}(j+1).Xs(2:end)-Tr.CVTwists{i}(j+1).Xs(1:end-1))*Lscale;
        xi_star_Z1 = reshape(Tr.CVTwists{i}(j+1).xi_star(:,2),6,nip);
        xi_star_Z2 = reshape(Tr.CVTwists{i}(j+1).xi_star(:,3),6,nip);
    
        for ii=1:nip-1
            h_a      = h(ii);
            Phi_Z1_a = Tr.PhiE_Z1(6*(isig-1)+1:6*isig,:);
            Phi_Z2_a = Tr.PhiE_Z2(6*(isig-1)+1:6*isig,:);
            xi_star_Z1_a = xi_star_Z1(:,ii); 
            xi_star_Z2_a = xi_star_Z2(:,ii); 
        
            [gstep((ij-1)*4+1:4*ij,:),S((ij-1)*6+1:6*ij,:),Sd((ij-1)*6+1:6*ij,:)]=JointQuantities_FD_mex(h_a,Phi_Z1_a,Phi_Z2_a,xi_star_Z1_a,xi_star_Z2_a,q,qd);
            Adgstepinv((ij-1)*6+1:6*ij,:)  = dinamico_Adjoint(ginv(gstep((ij-1)*4+1:4*ij,:)));
            ij = ij+1;
            isig = isig+1;
        end
        isig = isig+1;
    else
        isig = isig+1; %center of rigid link
    end
end
%% Forward Pass to compute g, eta, eta_plus, etadot, psi, etadot_plus, S_B, R, Q, Y, R_B, Q_B, Y_B
%initialize

g   = zeros(4*nsig,4);
eta = zeros(6*nsig,1); %J*qd
psi = zeros(6*nsig,1); %Jd*qd

S_B  = zeros(6*nsig,ndof); %Jacobian
Sd_B = zeros(6*nsig,ndof); %Derivative of Jacobian

g_tip   = repmat(eye(4),N,1);
S_Btip  = repmat(zeros(6,ndof),N,1);
Sd_Btip = repmat(zeros(6,ndof),N,1);
eta_tip = zeros(N*6,1);
psi_tip = zeros(N*6,1);

iLpre   = Tr.iLpre;
g_ini   = Tr.g_ini;

isig = 1;
ij  = 1;
for i=1:N
    if iLpre(i)>0 %imp for branched chain
        g_here       = g_tip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        S_Bhere      = Ad_g_ini_inv*S_Btip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        Sd_Bhere     = Ad_g_ini_inv*Sd_Btip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        psi_here     = Ad_g_ini_inv*psi_tip((iLpre(i)-1)*6+1:iLpre(i)*6);
        eta_here     = Ad_g_ini_inv*eta_tip((iLpre(i)-1)*6+1:iLpre(i)*6);
    else 
        g_here   = g_ini((i-1)*4+1:i*4,:);
        S_Bhere  = zeros(6,ndof);
        Sd_Bhere = zeros(6,ndof);
        psi_here = zeros(6,1);
        eta_here = zeros(6,1);
    end

    g((isig-1)*4+1:4*isig,:) = g_here;
    S_B((isig-1)*6+1:6*isig,:) = S_Bhere;
    Sd_B((isig-1)*6+1:6*isig,:) = Sd_Bhere;
    eta((isig-1)*6+1:6*isig) = eta_here;
    psi((isig-1)*6+1:6*isig) = psi_here; % to find FD

    % from 0 to 1 of rigid joint
    g_here = g_here*gstep((ij-1)*4+1:4*ij,:);
    S_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(S_Bhere+S((ij-1)*6+1:6*ij,:));
    Sd_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(Sd_Bhere+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,:)+Sd((ij-1)*6+1:6*ij,:));
    psi_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(psi_here+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,:)*qd+Sd((ij-1)*6+1:6*ij,:)*qd);
    eta_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(eta_here+S((ij-1)*6+1:6*ij,:)*qd);
    

    isig = isig+1;
    ij = ij+1;

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        %joint to CM
        gi        = Tr.VLinks(Tr.LinkIndex(i)).gi; 
        g_here    = g_here*gi;
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        S_Bhere   = Ad_gi_inv*S_Bhere;
        Sd_Bhere  = Ad_gi_inv*Sd_Bhere;
        psi_here  = Ad_gi_inv*psi_here;
        eta_here  = Ad_gi_inv*eta_here;

        g((isig-1)*4+1:4*isig,:) = g_here;
        S_B((isig-1)*6+1:6*isig,:) = S_Bhere;
        Sd_B((isig-1)*6+1:6*isig,:) = Sd_Bhere;
        psi((isig-1)*6+1:6*isig) = psi_here; % to find FD
        eta((isig-1)*6+1:6*isig) = eta_here;

        isig = isig+1;

        %CM to tip
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf; 
        g_here    = g_here*gf;
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        S_Bhere   = Ad_gf_inv*S_Bhere;
        Sd_Bhere  = Ad_gf_inv*Sd_Bhere;
        psi_here  = Ad_gf_inv*psi_here;
        eta_here  = Ad_gf_inv*eta_here;
        
    else
        %gi and gf of soft link are identity hence no transformation is done here
        %at X=0
        g((isig-1)*4+1:4*isig,:) = g_here;
        S_B((isig-1)*6+1:6*isig,:) = S_Bhere;
        Sd_B((isig-1)*6+1:6*isig,:) = Sd_Bhere;
        psi((isig-1)*6+1:6*isig) = psi_here; % to find FD
        eta((isig-1)*6+1:6*isig) = eta_here;
        
        isig = isig+1;
        
        nip = Tr.CVTwists{i}(j+1).nip;
        for ii = 1:nip-1
            g_here = g_here*gstep((ij-1)*4+1:4*ij,:);
            S_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(S_Bhere+S((ij-1)*6+1:6*ij,:));
            Sd_Bhere = Adgstepinv((ij-1)*6+1:6*ij,:)*(Sd_Bhere+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,:)+Sd((ij-1)*6+1:6*ij,:));
            psi_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(psi_here+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,:)*qd+Sd((ij-1)*6+1:6*ij,:)*qd);
            eta_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(eta_here+S((ij-1)*6+1:6*ij,:)*qd);
            
            ij = ij+1;

            g((isig-1)*4+1:4*isig,:) = g_here;
            S_B((isig-1)*6+1:6*isig,:) = S_Bhere;
            Sd_B((isig-1)*6+1:6*isig,:) = Sd_Bhere;
            psi((isig-1)*6+1:6*isig) = psi_here; % to find FD
            eta((isig-1)*6+1:6*isig) = eta_here;
            
            isig = isig+1;
        end

    end
    g_tip((i-1)*4+1:i*4,:)   = g_here;
    S_Btip((i-1)*6+1:i*6,:)  = S_Bhere;
    Sd_Btip((i-1)*6+1:i*6,:) = Sd_Bhere;
    psi_tip((i-1)*6+1:i*6) = psi_here;
    eta_tip((i-1)*6+1:i*6) = eta_here;
    

end

%% Backward Pass, compute the rest
%no need to save all. updated after current value is used. May have to save
%some for branched chain (later)

isig = nsig;
ij = nj;
ID_qqd = zeros(ndof,1);
M = zeros(ndof,ndof); %dID_dqdd

M_C_base = repmat(zeros(6,6),N,1);
F_C_qqd_base = repmat(zeros(6,1),N,1);
W_S_base = repmat(zeros(6,ndof),N,1);

for i=N:-1:1 %backwards

    M_C = zeros(6,6);
    F_C_qqd = zeros(6,1);
    W_S = zeros(6,ndof);
    for ic = Tr.iLnext{i} %fun if i has any children
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((ic-1)*4+1:ic*4,:)));
        coAd_g_ini = Ad_g_ini_inv';

        M_C     = M_C+coAd_g_ini*M_C_base(6*(ic-1)+1:6*ic,:)*Ad_g_ini_inv;
        F_C_qqd = F_C_qqd+coAd_g_ini*F_C_qqd_base(6*(ic-1)+1:6*ic,:);
        W_S     = W_S+coAd_g_ini*W_S_base(6*(ic-1)+1:6*ic,:);
    end

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        %tip to CM
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf; 
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        coAd_gf   = Ad_gf_inv';
        
        M_C     = coAd_gf*M_C*Ad_gf_inv;
        F_C_qqd = coAd_gf*F_C_qqd;
        W_S     = coAd_gf*W_S; %S_alpha is zero (imagine as a fixed joint)
        isig = isig-1;

        %CM to base
        gi        = Tr.VLinks(Tr.LinkIndex(i)).gi; 
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        coAd_gi   = Ad_gi_inv';
        
        M_ap1 = Tr.VLinks(Tr.LinkIndex(i)).M; %at alpha + 1. 
        F_ap1_qqd = -M_ap1*dinamico_Adjoint(ginv(g(isig*4+1:4*(isig+1),:)))*Tr.G+M_ap1*psi(isig*6+1:6*(isig+1))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1*eta(isig*6+1:6*(isig+1));

        M_C     = coAd_gi*(M_ap1+M_C)*Ad_gi_inv;
        F_C_qqd = coAd_gi*(F_ap1_qqd+F_C_qqd);
        W_S     = coAd_gi*W_S; %S_alpha is zero (imagine as a fixed joint)

        %for rigid joint
        M_ap1 = zeros(6,6); %at base
        F_ap1_qqd = zeros(6,1);

    else
        %assuming gf and gi of soft are identities
        Ws = Tr.CVTwists{i}(j+1).Ws; %L*w_gauss
        nip = Tr.CVTwists{i}(j+1).nip;
        isig = isig-1;
        for ii=nip-1:-1:1
        
            coAdgstep = Adgstepinv((ij-1)*6+1:6*ij,:)';
        
            M_ap1 = Ws(ii+1)*Tr.CVTwists{i}(2).Ms(ii*6+1:6*(ii+1),:); %at alpha + 1. multiplied with weight
            F_ap1_qqd = -M_ap1*dinamico_Adjoint(ginv(g(isig*4+1:4*(isig+1),:)))*Tr.G+M_ap1*psi(isig*6+1:6*(isig+1))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1*eta(isig*6+1:6*(isig+1));
            
            M_C = coAdgstep*(M_ap1+M_C)*Adgstepinv((ij-1)*6+1:6*ij,:);
            F_C_qqd = coAdgstep*(F_ap1_qqd+F_C_qqd);
            W_S = M_C*S((ij-1)*6+1:6*ij,:)+coAdgstep*W_S;
        
            ID_qqd  = ID_qqd+S((ij-1)*6+1:6*ij,:)'*F_C_qqd;
            M       = M+S((ij-1)*6+1:6*ij,:)'*(M_C*S_B((isig-1)*6+1:6*isig,:)+W_S);

            isig = isig-1;
            ij = ij-1;
        
        end
        %for rigid joint
        M_ap1 = Ws(1)*Tr.CVTwists{i}(2).Ms(1:6,:); %at X = 0 usually Ws is 0
        F_ap1_qqd = -M_ap1*dinamico_Adjoint(ginv(g(isig*4+1:4*(isig+1),:)))*Tr.G+M_ap1*psi(isig*6+1:6*(isig+1))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1*eta(isig*6+1:6*(isig+1));
    end
    %joint
    coAdgstep = Adgstepinv((ij-1)*6+1:6*ij,:)';
    
    M_C = coAdgstep*(M_ap1+M_C)*Adgstepinv((ij-1)*6+1:6*ij,:);
    F_C_qqd = coAdgstep*(F_ap1_qqd+F_C_qqd);
    W_S = M_C*S((ij-1)*6+1:6*ij,:)+coAdgstep*W_S;

    ID_qqd  = ID_qqd+S((ij-1)*6+1:6*ij,:)'*F_C_qqd;
    M       = M+S((ij-1)*6+1:6*ij,:)'*(M_C*S_B((isig-1)*6+1:6*isig,:)+W_S);
    
    isig = isig-1;
    ij = ij-1;

    M_C_base((i-1)*6+1:i*6,:) = M_C;
    F_C_qqd_base((i-1)*6+1:i*6) = F_C_qqd;
    W_S_base((i-1)*6+1:i*6,:) = W_S;

end

%% derivatives of internal wrenches
tau = -Tr.K*q-Tr.D*qd; %K and D are same even if we use extended basis, modify only for POD case
B = Tr.Bqj1; %only 1 DOF joints
%% Combine

%find M_new B_new
u = zeros(Tr.nact,1);

for i=1:Tr.n_jact
    if ~Tr.WrenchControlled(i)  %swapping for joint cooridnate controlled joints
        M_temp          = M;
        M(:,Tr.i_jactq(i)) = -B(:,i); 
        B(:,i)          = -M_temp(:,Tr.i_jactq(i));  %Tau=Bq*u
        u(i)             = qqdqddjt(i,3);
    end
end
tau = tau+B*u;

qdd = M\(tau-ID_qqd);
u_act = zeros(Tr.nact,1);
for i=1:Tr.n_jact
    if ~Tr.WrenchControlled(i)
        u_act(i) = qdd(Tr.i_jactq(i));
        qdd(Tr.i_jactq(i)) = u(i); %replace u with qdd if wrench controlled
    end
end
qddu = [qdd;u_act];
end