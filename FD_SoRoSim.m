%Function to calculate FD for a given q, qd, and u using ABA, SoRoSim
%Algorithm
function qdd=FD_SoRoSim(Tr,t,qqd,uqt)
%disp(['Current time t: ', num2str(t)]);

ndof   = Tr.ndof;
q   = qqd(1:Tr.ndof);
qd  = qqd(Tr.ndof+1:2*Tr.ndof);

N      = Tr.N;

qqdqddjt = uqt(t);

for i=1:Tr.n_jact
    if ~Tr.WrenchControlled(i)
        q(Tr.i_jactq(i))  = qqdqddjt(i,1); %replacing q_k and qdot_k
        qd(Tr.i_jactq(i)) = qqdqddjt(i,2);
    end
end

M = zeros(ndof,ndof);
F = zeros(ndof,1);

G = Tr.G; %gravity

dof_start = 1; %starting dof of current piece

g_ini     = Tr.g_ini; %initial configuration of all link wrt its previous link
g_Ltip    = repmat(eye(4),N,1);
J_Ltip    = repmat(zeros(6,ndof),N,1);
eta_Ltip  = zeros(N*6,1); %total velocity J*qd+eta_t
psi_Ltip  = zeros(N*6,1);

iLpre     = Tr.iLpre;

for i=1:N

    if iLpre(i)>0
        g_here       = g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
        Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
        J_here       = Ad_g_ini_inv*J_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        
        eta_here   = Ad_g_ini_inv*eta_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6);
        psi_here   = Ad_g_ini_inv*psi_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6);
    else
        g_here   = g_ini((i-1)*4+1:i*4,:);
        J_here   = zeros(6,ndof);
        
        eta_here   = zeros(6,1);
        psi_here   = zeros(6,1);
    end
    
    %Joint
    dof_here = Tr.CVTwists{i}(1).dof;
    q_here   = q(dof_start:dof_start+dof_here-1);
    qd_here  = qd(dof_start:dof_start+dof_here-1);
    B_here   = Tr.CVTwists{i}(1).B;
    xi_star  = Tr.CVTwists{i}(1).xi_star;

    if dof_here==0 %fixed joint (N)
        g_joint   = eye(4);
        TgB_here  = zeros(6,ndof);
        TgBd_here = zeros(6,ndof);
    else
        xi          = B_here*q_here+xi_star;
        xid         = B_here*qd_here;

        [g_joint,Tg,Tgd]=variable_expmap_gTgTgd_mex(xi,xid);

        TgB_here                                    = zeros(6,ndof);
        TgB_here(:,dof_start:dof_start+dof_here-1)  = Tg*B_here;
        TgBd_here                                   = zeros(6,ndof);
        TgBd_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*Tg*B_here+Tgd*B_here;
    end

    %updating g, Jacobian, Jacobian_dot and eta
    g_here         = g_here*g_joint;
    Ad_g_joint_inv = dinamico_Adjoint(ginv(g_joint));
    J_here         = Ad_g_joint_inv*(J_here+TgB_here);
    eta_here       = Ad_g_joint_inv*(eta_here+TgB_here(:,dof_start:dof_start+dof_here-1)*qd_here);
    psi_here       = Ad_g_joint_inv*(psi_here+TgBd_here(:,dof_start:dof_start+dof_here-1)*qd_here);


    if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'

        gi        = Tr.VLinks(Tr.LinkIndex(i)).gi;
        g_here    = g_here*gi;
        Ad_gi_inv = dinamico_Adjoint(ginv(gi));
        J_here    = Ad_gi_inv*J_here;
        eta_here  = Ad_gi_inv*eta_here;
        psi_here  = Ad_gi_inv*psi_here;

        
        M_here = Tr.VLinks(Tr.LinkIndex(i)).M;

        F = F+J_here'*(M_here*dinamico_Adjoint(ginv(g_here))*G-M_here*psi_here-dinamico_coadj(eta_here)*M_here*eta_here);
        M = M+J_here'*M_here*J_here;
        
        % bringing all quantities to the end of rigid link
        gf        = Tr.VLinks(Tr.LinkIndex(i)).gf;
        g_here    = g_here*gf;
        Ad_gf_inv = dinamico_Adjoint(ginv(gf));
        J_here    = Ad_gf_inv*J_here;
        eta_here  = Ad_gf_inv*eta_here;
        psi_here  = Ad_gf_inv*psi_here;
        
    end

    dof_start = dof_start+dof_here;

    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1 %will run only if soft link

        dof_here   = Tr.CVTwists{i}(j+1).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        qd_here    = qd(dof_start:dof_start+dof_here-1);
        
        Lscale  = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
        
        xi_star = Tr.CVTwists{i}(j+1).xi_star;
        Ms      = Tr.CVTwists{i}(j+1).Ms;
        Xs      = Tr.CVTwists{i}(j+1).Xs;
        Ws      = Tr.CVTwists{i}(j+1).Ws;
        nip     = Tr.CVTwists{i}(j+1).nip;
       
        ii = 1;
        if Ws(ii)>0
            M_here = Ws(ii)*Ms(6*(ii-1)+1:6*ii,:);

            F = F+J_here'*(M_here*dinamico_Adjoint(ginv(g_here))*G-M_here*psi_here-dinamico_coadj(eta_here)*M_here*eta_here);
            M = M+J_here'*M_here*J_here; 

        end

        for ii=2:nip

            H = (Xs(ii)-Xs(ii-1))*Lscale;
                
            xi_Z1here = xi_star(6*(ii-2)+1:6*(ii-1),2); 
            xi_Z2here = xi_star(6*(ii-2)+1:6*(ii-1),3);
                
            B_Z1here  = Tr.CVTwists{i}(j+1).B_Z1(6*(ii-2)+1:6*(ii-1),:);%note this step
            B_Z2here  = Tr.CVTwists{i}(j+1).B_Z2(6*(ii-2)+1:6*(ii-1),:);

            xi_Z1here = B_Z1here*q_here+xi_Z1here;
            xi_Z2here = B_Z2here*q_here+xi_Z2here;

            xid_Z1here  = B_Z1here*qd_here;
            ad_xi_Z1here = dinamico_adj(xi_Z1here);

            BGamma_here  = (H/2)*(B_Z1here+B_Z2here)+... %dBqdq = B
                               ((sqrt(3)*H^2)/12)*(ad_xi_Z1here*B_Z2here-dinamico_adj(xi_Z2here)*B_Z1here);

            Gammadd_Z4_dq_here = ((sqrt(3)*H^2)/6)*dinamico_adj(xid_Z1here)*B_Z2here; 
            Gammad_here   = BGamma_here*qd_here;

            Gamma_here = (H/2)*(xi_Z1here+xi_Z2here)+...
                         ((sqrt(3)*H^2)/12)*ad_xi_Z1here*xi_Z2here;  


            [gh,TGamma_here,TGammad_here] = variable_expmap_gTgTgd_mex(Gamma_here,Gammad_here); % mex code, C program

            TBGamma_here                                    = zeros(6,ndof);
            TBGamma_here(:,dof_start:dof_start+dof_here-1)  = TGamma_here*BGamma_here;
            TBGammad_here                                   = zeros(6,ndof);
            TBGammad_here(:,dof_start:dof_start+dof_here-1) = dinamico_adj(eta_here)*TBGamma_here(:,dof_start:dof_start+dof_here-1)+TGammad_here*BGamma_here;
            TBGammad_here(:,dof_start:dof_start+dof_here-1) = TBGammad_here(:,dof_start:dof_start+dof_here-1)+TGamma_here*Gammadd_Z4_dq_here;
            
            %updating g, Jacobian, Jacobian_dot and eta
            g_here     = g_here*gh;
            Ad_gh_inv  = dinamico_Adjoint(ginv(gh));
            J_here     = Ad_gh_inv*(J_here+TBGamma_here); %full
            eta_here   = Ad_gh_inv*(eta_here+TBGamma_here(:,dof_start:dof_start+dof_here-1)*qd_here);
            psi_here   = Ad_gh_inv*(psi_here+TBGammad_here(:,dof_start:dof_start+dof_here-1)*qd_here);
            
            
            %integrals evaluation
            if Ws(ii)>0
                M_here = Ws(ii)*Ms(6*(ii-1)+1:6*ii,:);

                F = F+J_here'*(M_here*dinamico_Adjoint(ginv(g_here))*G-M_here*psi_here-dinamico_coadj(eta_here)*M_here*eta_here); %rescale
                M = M+J_here'*M_here*J_here;  
            end
        end

        dof_start = dof_start+dof_here;
    end
    
    g_Ltip((i-1)*4+1:i*4,:) = g_here;
    J_Ltip((i-1)*6+1:i*6,:) = J_here;
    eta_Ltip((i-1)*6+1:i*6) = eta_here;
    psi_Ltip((i-1)*6+1:i*6) = psi_here;

end
%% derivatives of internal wrenches
tau = -Tr.K*q-Tr.D*qd; %K and D are same even if we use extended basis, modify only for POD case
B = Tr.Bqj1; %only 1 DOF joints
%% Combine

%find M_new B_new
u    = zeros(Tr.nact,1);
for i=1:Tr.n_jact
    if ~Tr.WrenchControlled(i)  %swapping for joint cooridnate controlled joints
        M_temp          = M;
        M(:,Tr.i_jactq(i)) = -B(:,i); 
        B(:,i)          = -M_temp(:,Tr.i_jactq(i));  %Tau=Bq*u
        u(i)            = qqdqddjt(i,3);
    end
end
tau = tau+B*u;

qdd = M\(tau+F);
for i=1:Tr.n_jact
    if ~Tr.WrenchControlled(i)
        qdd(Tr.i_jactq(i)) = u(i); %replace u with qdd if wrench controlled
    end
end

end