function Res = DAE(Tr,t,x,xd,uqt)
N    = Tr.N;
ndof = Tr.ndof;
nsig = Tr.nsig; %nip+1 for soft link (rigid joint+soft body) and 2 for rigid link
nj   = nsig-N; %nip joints for soft link (rigid joint+distributed joint) and 1 for rigid link

q   = x(1:Tr.ndof);
qd  = x(Tr.ndof+1:2*Tr.ndof);
u   = x(2*Tr.ndof+1:end);

qd_xd = xd(1:Tr.ndof);
qdd = xd(Tr.ndof+1:2*Tr.ndof);

%% Evaluate Quantities at every joints

h = zeros(1,nj);
gstep = zeros(4*nj,4);
Adgstepinv = zeros(6*nj,6);
S = zeros(6*nj,ndof);
Sd = zeros(6*nj,ndof);


j = 1; % only 1 soft division
ij = 1;
isig = 1;

for i=1:N
    h(ij) = 1;
    Phi = Tr.PhiE(6*(isig-1)+1:6*isig,:);
    [gstep((ij-1)*4+1:4*ij,:),S((ij-1)*6+1:6*ij,:),Sd((ij-1)*6+1:6*ij,:)]=JointQuantitiesR_FD_mex(Phi,q,qd);
    Adgstepinv((ij-1)*6+1:6*ij,:)  = dinamico_Adjoint(ginv(gstep((ij-1)*4+1:4*ij,:)));
    ij = ij+1;
    isig = isig+1;
    %for soft
    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='s'

        nip = Tr.CVTwists{i}(j+1).nip;
        Lscale = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        h(ij:ij+nip-2) = (Tr.CVTwists{i}(j+1).Xs(2:end)-Tr.CVTwists{i}(j+1).Xs(1:end-1))*Lscale;
        xi_star_Z1 = reshape(Tr.CVTwists{i}(j+1).xi_star(:,2),6,nip);
        xi_star_Z2 = reshape(Tr.CVTwists{i}(j+1).xi_star(:,3),6,nip);
    
        for ii=1:nip-1
            Phi_Z1_ij = Tr.PhiE_Z1(6*(isig-1)+1:6*isig,:);
            Phi_Z2_ij = Tr.PhiE_Z2(6*(isig-1)+1:6*isig,:);
            xi_star_Z1_ij = xi_star_Z1(:,ii); 
            xi_star_Z2_ij = xi_star_Z2(:,ii); 

            [gstep((ij-1)*4+1:4*ij,:),S((ij-1)*6+1:6*ij,:),Sd((ij-1)*6+1:6*ij,:)]=JointQuantities_FD_mex(h(ij),Phi_Z1_ij,Phi_Z2_ij,xi_star_Z1_ij,xi_star_Z2_ij,q,qd);
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
etadot = zeros(6*nsig,1); %Jd*qd

g_tip   = repmat(eye(4),N,1);

eta_tip = zeros(N*6,1);
etadot_tip = zeros(N*6,1);

iLpre   = Tr.iLpre;
g_ini   = Tr.g_ini;

isig = 1;
ij  = 1;
for i=1:N
    if iLpre(i)>0 %imp for branched chain
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));

        g_here       = g_tip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        etadot_here  = Ad_g_ini_inv*etadot_tip((iLpre(i)-1)*6+1:iLpre(i)*6);
        eta_here     = Ad_g_ini_inv*eta_tip((iLpre(i)-1)*6+1:iLpre(i)*6);
    else 

        g_here   = g_ini((i-1)*4+1:i*4,:);
        etadot_here = zeros(6,1);
        eta_here    = zeros(6,1);
    end
    
    %0 of joint

    g((isig-1)*4+1:4*isig,:) = g_here;
    etadot((isig-1)*6+1:6*isig) = etadot_here;
    eta((isig-1)*6+1:6*isig) = eta_here;

    eta_plus_here    = eta_here+S((ij-1)*6+1:6*ij,:)*qd;
    etadot_plus_here = etadot_here+S((ij-1)*6+1:6*ij,:)*qdd+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,:)*qd+Sd((ij-1)*6+1:6*ij,:)*qd;
    
    % from 0 to 1 of rigid joint

    g_here = g_here*gstep((ij-1)*4+1:4*ij,:);
    etadot_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(etadot_plus_here);
    eta_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(eta_plus_here);

    isig = isig+1;
    ij = ij+1;

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        %joint to CM
        gi        = Tr.VLinks(Tr.LinkIndex(i)).gi; 
        
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));

        g_here    = g_here*gi;
        etadot_here  = Ad_gi_inv*etadot_here;
        eta_here  = Ad_gi_inv*eta_here;
        
        g((isig-1)*4+1:4*isig,:) = g_here;
        etadot((isig-1)*6+1:6*isig) = etadot_here;
        eta((isig-1)*6+1:6*isig) = eta_here; % to find FD

        isig = isig+1;

        %CM to tip
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf; 

        Ad_gf_inv = dinamico_Adjoint(ginv(gf));

        g_here    = g_here*gf;
        etadot_here  = Ad_gf_inv*etadot_here;
        eta_here  = Ad_gf_inv*eta_here;
        
    else
        %gi and gf of soft link are identity hence no transformation is done here
        %at X=0
        
        g((isig-1)*4+1:4*isig,:) = g_here;
        etadot((isig-1)*6+1:6*isig) = etadot_here;
        eta((isig-1)*6+1:6*isig) = eta_here; % to find FD
        
        isig = isig+1;
        
        nip = Tr.CVTwists{i}(j+1).nip;
        for ii = 1:nip-1

            eta_plus_here    = eta_here+S((ij-1)*6+1:6*ij,:)*qd;
            etadot_plus_here = etadot_here+S((ij-1)*6+1:6*ij,:)*qdd+dinamico_adj(eta_here)*S((ij-1)*6+1:6*ij,:)*qd+Sd((ij-1)*6+1:6*ij,:)*qd;

            g_here = g_here*gstep((ij-1)*4+1:4*ij,:);
            etadot_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(etadot_plus_here);
            eta_here = Adgstepinv((ij-1)*6+1:6*ij,:)*(eta_plus_here);
    
            g((isig-1)*4+1:4*isig,:) = g_here;
            etadot((isig-1)*6+1:6*isig) = etadot_here;
            eta((isig-1)*6+1:6*isig) = eta_here; % to find FD 
            
            isig = isig+1;
            ij = ij+1;
        end

    end
    % saving tip values

    g_tip((i-1)*4+1:i*4,:)   = g_here;
    etadot_tip((i-1)*6+1:i*6) = etadot_here;
    eta_tip((i-1)*6+1:i*6) = eta_here;
    
end

%% Backward Pass, compute the rest
%no need to save all. updated after current value is used.

isig = nsig;
ij = nj;

ID = zeros(ndof,1);

F_C_base = repmat(zeros(6,1),N,1);

for i=N:-1:1 %backwards

    F_C = zeros(6,1);

    for ic = Tr.iLnext{i} %runs if i has any children. iLnext can be avoided respecting sorosim's order of linkage creation: project all to gf_parent before
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((ic-1)*4+1:ic*4,:)));
        coAd_g_ini = Ad_g_ini_inv';

        F_C = F_C+coAd_g_ini*F_C_base(6*(ic-1)+1:6*ic);
    end

    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
        %tip to CM
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf; 
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        coAd_gf   = Ad_gf_inv';
        
        F_C = coAd_gf*F_C;
        isig = isig-1;

        %CM to base
        gi        = Tr.VLinks(Tr.LinkIndex(i)).gi; 
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        coAd_gi   = Ad_gi_inv';
        
        M_ap1 = Tr.VLinks(Tr.LinkIndex(i)).M; %at alpha + 1.
        F_ap1 = -M_ap1*dinamico_Adjoint(ginv(g(isig*4+1:4*(isig+1),:)))*Tr.G+M_ap1*etadot(isig*6+1:6*(isig+1))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1*eta(isig*6+1:6*(isig+1));

        F_C = coAd_gi*(F_ap1+F_C);

        %for rigid joint
        F_ap1 = zeros(6,1);
    else
        %assuming gf and gi of soft are identities
        Ws = Tr.CVTwists{i}(j+1).Ws; %L*w_gauss
        nip = Tr.CVTwists{i}(j+1).nip;
        isig = isig-1;
        for ii=nip-1:-1:1
        
            coAdgstep = Adgstepinv((ij-1)*6+1:6*ij,:)';
        
            M_ap1 = Ws(ii+1)*Tr.CVTwists{i}(2).Ms(ii*6+1:6*(ii+1),:); %at alpha + 1. multiplied with weight
            F_ap1 = -M_ap1*dinamico_Adjoint(ginv(g(isig*4+1:4*(isig+1),:)))*Tr.G+M_ap1*etadot(isig*6+1:6*(isig+1))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1*eta(isig*6+1:6*(isig+1));

            F_C = coAdgstep*(F_ap1+F_C);
            
            ID = ID + S((ij-1)*6+1:6*ij,:)'*F_C;

            isig = isig-1;
            ij = ij-1;
        
        end
        %for rigid joint
        M_ap1 = Ws(1)*Tr.CVTwists{i}(2).Ms(1:6,:); %at X = 0 usually Ws is 0
        F_ap1 = -M_ap1*dinamico_Adjoint(ginv(g(isig*4+1:4*(isig+1),:)))*Tr.G+M_ap1*etadot(isig*6+1:6*(isig+1))+dinamico_coadj(eta(isig*6+1:6*(isig+1)))*M_ap1*eta(isig*6+1:6*(isig+1));
    end
    %joint
    coAdgstep = Adgstepinv((ij-1)*6+1:6*ij,:)';
    
    F_C = coAdgstep*(F_ap1+F_C);

    ID = ID + S((ij-1)*6+1:6*ij,:)'*F_C;
    
    isig = isig-1;
    ij = ij-1;

    F_C_base((i-1)*6+1:i*6) = F_C;

end

%% derivatives of internal wrenches

B = Tr.Bqj1; %only 1 DOF joints, its derivative is 0
tau = -Tr.K*q-Tr.D*qd+B*u;

qqdqddjt = uqt(t);
gc = B'*q-qqdqddjt(:,1); %replace u with qdd if wrench controlled

Res = [qd-qd_xd;tau-ID;gc];
end