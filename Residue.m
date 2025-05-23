function Res = Residue(Tr,q,u) %Res =  ID(q) - tau(q,u)

ndof = Tr.ndof;

i=1;
j=1;
nip = Tr.CVTwists{i}(j+1).nip;
Lscale = Tr.VLinks(Tr.LinkIndex(i)).lp{j};

%% Evaluate Quantities at every joints

h = (Tr.CVTwists{i}(j+1).Xs(2:end)-Tr.CVTwists{i}(j+1).Xs(1:end-1))*Lscale;

xi_star_Z1 = reshape(Tr.CVTwists{i}(j+1).xi_star(:,2),6,nip);
xi_star_Z2 = reshape(Tr.CVTwists{i}(j+1).xi_star(:,3),6,nip);


gstep = zeros(4*(nip-1),4);
Adgstepinv = zeros(6*(nip-1),6);
S = zeros(6*(nip-1),ndof);

for a=1:nip-1
    h_a      = h(a);
    Phi_Z1_a = Tr.CVTwists{1}(2).B_Z1((a-1)*6+1:6*a,:);
    Phi_Z2_a = Tr.CVTwists{1}(2).B_Z2((a-1)*6+1:6*a,:);
    xi_star_Z1_a = xi_star_Z1(:,a); 
    xi_star_Z2_a = xi_star_Z2(:,a); 

    [~,~,gstep((a-1)*4+1:4*a,:),~,S((a-1)*6+1:6*a,:),~,~,~]=JointQuantities_statics_mex(h_a,Phi_Z1_a,Phi_Z2_a,xi_star_Z1_a,xi_star_Z2_a,q);
    Adgstepinv((a-1)*6+1:6*a,:)  = dinamico_Adjoint(ginv(gstep((a-1)*4+1:4*a,:)));
end

%% Forward Pass to compute g, eta, eta_plus, etadot, psi, etadot_plus, S_B, R, Q, Y, R_B, Q_B, Y_B
%initialize

g   = zeros(4*nip,4);

a = 1;
g((a-1)*4+1:4*a,:) = Tr.g_ini((i-1)*4+1:i*4,:);

for a = 1:nip-1
    % find quantities at 2 to nip (for g, eta, and etadot and S_B, R_B, Q_B, and Y_B)
    g(a*4+1:4*(a+1),:) = g((a-1)*4+1:4*a,:)*gstep((a-1)*4+1:4*a,:);
end

%% Backward Pass, compute the rest
%no need to save all. updated after current value is used. May have to save
%some for branched chain (later)

F_C = zeros(6,1);
ID = zeros(ndof,1);

Ws = Tr.CVTwists{i}(j+1).Ws; %L*w_gauss
for a=nip-1:-1:1
    coAdgstep = Adgstepinv((a-1)*6+1:6*a,:)';   
    M_ap1 = Ws(a+1)*Tr.CVTwists{1}(2).Ms(a*6+1:6*(a+1),:); %at alpha + 1. multiplied with weight
    F_ap1 = -M_ap1*dinamico_Adjoint(ginv(g(a*4+1:4*(a+1),:)))*Tr.G;
    F_C = coAdgstep*(F_ap1+F_C);
    ID     = ID+S((a-1)*6+1:6*a,:)'*(F_C);
end

%% derivatives of internal wrenches
tau = -Tr.K*q;

for ii=1:nip %could be combined with parallel processing or one of previous loop
    if Ws(ii)>0
        Phi = Tr.CVTwists{1}(2).B((ii-1)*6+1:ii*6,:);
        xi  = Phi*q+Tr.CVTwists{i}(j+1).xi_star((ii-1)*6+1:ii*6,1);
        Phi_au = zeros(6,1);
        xihat_123  = [0 -xi(3) xi(2) xi(4);xi(3) 0 -xi(1) xi(5);-xi(2) xi(1) 0 xi(6)];%4th row is avoided to speedup calculation
        for ia=1:Tr.n_sact
            dc  = Tr.dc{ia,i}{j}(:,ii);
            dcp = Tr.dcp{ia,i}{j}(:,ii);
            Phi_au_ia = SoftActuator_FD_mex(u(ia),dc,dcp,xihat_123);
            Phi_au = Phi_au+Phi_au_ia;
        end
        tau = tau+Ws(ii)*Phi'*Phi_au;
    end
end

Res = ID-tau;
end