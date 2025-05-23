%Function for the plot of dynamic simulation
%Last modified by Anup Teejo Mathew - 25/05/2021
% Edited by Yuki - 15/04/2025
function plotqqd_test(Tr,t,qqd,i,xyz_gt_pos,xyz_pred_pos,out_folder)
close all

PlottingParameters = Tr.PlotParameters;

N         = Tr.N;
g_ini     = Tr.g_ini;
iLpre     = Tr.iLpre;

tic
tmax        = max(t);
% 保存视频到指定输出文件夹
video_filename = fullfile(out_folder, ['predict_', num2str(i), '.avi']);
v = VideoWriter(video_filename);
FrameRate = PlottingParameters.FrameRateValue * 2;
v.FrameRate = FrameRate;
open(v);


hfig      = figure('units','normalized','outerposition',[0 0 1 1]);
hfig.WindowState = 'maximized';
set(gca,'CameraPosition',PlottingParameters.CameraPosition,...
        'CameraTarget',PlottingParameters.CameraTarget,...
        'CameraUpVector',PlottingParameters.CameraUpVector,...
        'FontSize',24);
    
if PlottingParameters.Light
    camlight(PlottingParameters.Az_light,PlottingParameters.El_light)
end
%     view(45,45);
axis equal
grid on
hold on
xlabel('x (m)')
ylabel('y (m)')
zlabel('z (m)') 
axis([PlottingParameters.X_lim PlottingParameters.Y_lim PlottingParameters.Z_lim]);

set(get(gca, 'Title'), 'Interpreter', 'latex');
set(get(gca, 'XLabel'), 'Interpreter', 'latex');
set(get(gca, 'YLabel'), 'Interpreter', 'latex');
set(get(gca, 'ZLabel'), 'Interpreter', 'latex');
set(gca, 'TickLabelInterpreter', 'latex');

% Adjust legend and text annotations to use LaTeX interpreter
hLegend = findobj(gcf, 'Type', 'Legend');
set(hLegend, 'Interpreter', 'latex');

hText = findobj(gca, 'Type', 'Text');
set(hText, 'Interpreter', 'latex');
set(gca,'FontSize',24)
%axis tight

drawnow
drawnow

lt = length(0:1/FrameRate:tmax);
ptip = NaN(3,lt);
it=1;

for tt=0:1/FrameRate:tmax

    delete(findobj('type', 'patch'));
    title(strcat('t= ',num2str(tt)))
   
    qqdtt = interp1(t,qqd,tt);
    q     = qqdtt(1:Tr.ndof)';
    qd     = qqdtt(Tr.ndof+1:2*Tr.ndof)';
    
    dof_start = 1;
    g_Ltip    = repmat(eye(4),N,1);
    J_Ltip    = repmat(zeros(6,Tr.ndof),N,1);
    
    for i=1:N % number of links

        if tt>0
            if exist('hv1','var')
            if isgraphics(hv1)
                delete(hv1);
            end
            end

            delete(h2);
        end
        
        if iLpre(i)>0
            g_here=g_Ltip((iLpre(i)-1)*4+1:iLpre(i)*4,:)*g_ini((i-1)*4+1:i*4,:);
            Ad_g_ini_inv = dinamico_Adjoint(ginv(g_ini((i-1)*4+1:i*4,:)));
            J_here       = Ad_g_ini_inv*J_Ltip((iLpre(i)-1)*6+1:iLpre(i)*6,:);
        else
            g_here=g_ini((i-1)*4+1:i*4,:);
            J_here   = zeros(6,Tr.ndof);
        end
        
        %joint
        dof_here   = Tr.CVTwists{i}(1).dof;
        q_here     = q(dof_start:dof_start+dof_here-1);
        B_here     = Tr.CVTwists{i}(1).B;
        xi_star    = Tr.CVTwists{i}(1).xi_star;

        if dof_here==0 %fixed joint (N)
            g_joint    = eye(4);
            TgB_here  = zeros(6,Tr.ndof);
        else
            xi         = B_here*q_here+xi_star;

            [g_joint,Tg]=variable_expmap_gTg_mex(xi);
            TgB_here                                    = zeros(6,Tr.ndof);
            TgB_here(:,dof_start:dof_start+dof_here-1)  = Tg*B_here;
        end
        g_here     = g_here*g_joint;
        Ad_g_joint_inv = dinamico_Adjoint(ginv(g_joint));
        J_here         = Ad_g_joint_inv*(J_here+TgB_here);
        
        n_r   = Tr.VLinks(Tr.LinkIndex(i)).n_r;
        if Tr.VLinks(Tr.LinkIndex(i)).CS=='R'
            n_r=5;
        end
        n_l   = Tr.VLinks(Tr.LinkIndex(i)).n_l;
        color = Tr.VLinks(Tr.LinkIndex(i)).color;
        
        if Tr.VLinks(Tr.LinkIndex(i)).L>0
        if Tr.VLinks(Tr.LinkIndex(i)).linktype=='r'
            L          = Tr.VLinks(Tr.LinkIndex(i)).L;
            gi         = Tr.VLinks(Tr.LinkIndex(i)).gi;
            g_here     = g_here*gi;
            Ad_gi_inv = dinamico_Adjoint(ginv(gi));
            J_here    = Ad_gi_inv*J_here;
            
            if ~Tr.VLinks(Tr.LinkIndex(i)).CPF
                Xr         = linspace(0,L,n_l);
                g_hereR    = g_here*[eye(3) [-Tr.VLinks(Tr.LinkIndex(i)).L/2;0;0];0 0 0 1]; 
                dx         = Xr(2)-Xr(1);

                Xpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
                Ypatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
                Zpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
                i_patch = 1;

                if Tr.VLinks(Tr.LinkIndex(i)).CS=='C'

                    r_fn  = Tr.VLinks(Tr.LinkIndex(i)).r;
                    r     = r_fn(0);
                    theta = linspace(0,2*pi,n_r);
                    x     = zeros(1,n_r);
                    y     = r*sin(theta);
                    z     = r*cos(theta);
                    pos   = [x;y;z;ones(1,n_r)];

                elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='R'

                    h_fn  = Tr.VLinks(Tr.LinkIndex(i)).h;
                    w_fn  = Tr.VLinks(Tr.LinkIndex(i)).w;
                    h     = h_fn(0);
                    w     = w_fn(0);
                    x     = [0 0 0 0 0];
                    y     = [h/2 -h/2 -h/2 h/2 h/2];
                    z     = [w/2 w/2 -w/2 -w/2 w/2];
                    pos   = [x;y;z;ones(1,5)];

                elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='E'

                    a_fn  = Tr.VLinks(Tr.LinkIndex(i)).a;
                    b_fn  = Tr.VLinks(Tr.LinkIndex(i)).b;
                    a     = a_fn(0);
                    b     = b_fn(0);
                    theta = linspace(0,2*pi,n_r);
                    x     = zeros(1,n_r);
                    y     = a*sin(theta);
                    z     = b*cos(theta);
                    pos   = [x;y;z;ones(1,n_r)];
                end

                pos_here = g_hereR*pos;
                x_here   = pos_here(1,:);
                y_here   = pos_here(2,:);
                z_here   = pos_here(3,:);

                Xpatch(:,i_patch) = x_here';
                Ypatch(:,i_patch) = y_here';
                Zpatch(:,i_patch) = z_here';
                i_patch           = i_patch+1;

                x_pre    = x_here;
                y_pre    = y_here;
                z_pre    = z_here;

                for ii=2:n_l

                    if Tr.VLinks(Tr.LinkIndex(i)).CS=='C'

                        r     = r_fn(Xr(ii)/L);
                        theta = linspace(0,2*pi,n_r);
                        x     = zeros(1,n_r);
                        y     = r*sin(theta);
                        z     = r*cos(theta);
                        pos   = [x;y;z;ones(1,n_r)];

                    elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='R'

                        h     = h_fn(Xr(ii)/L);
                        w     = w_fn(Xr(ii)/L);
                        x     = [0 0 0 0 0];
                        y     = [h/2 -h/2 -h/2 h/2 h/2];
                        z     = [w/2 w/2 -w/2 -w/2 w/2];
                        pos   = [x;y;z;ones(1,5)];

                    elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='E'

                        a     = a_fn(Xr(ii)/L);
                        b     = b_fn(Xr(ii)/L);
                        theta = linspace(0,2*pi,n_r);
                        x     = zeros(1,n_r);
                        y     = a*sin(theta);
                        z     = b*cos(theta);
                        pos   = [x;y;z;ones(1,n_r)];
                    end

                    g_hereR  = g_hereR*[eye(3) [dx;0;0];0 0 0 1];
                    pos_here = g_hereR*pos;
                    x_here   = pos_here(1,:);
                    y_here   = pos_here(2,:);
                    z_here   = pos_here(3,:);

                    %Plotting rigid link
                    for jj=1:n_r-1

                        Xpatch(1:5,i_patch)   = [x_pre(jj) x_here(jj) x_here(jj+1) x_pre(jj+1) x_pre(jj)]';
                        Xpatch(6:end,i_patch) = x_pre(jj)*ones(n_r-5,1);
                        Ypatch(1:5,i_patch)   = [y_pre(jj) y_here(jj) y_here(jj+1) y_pre(jj+1) y_pre(jj)]';
                        Ypatch(6:end,i_patch) = y_pre(jj)*ones(n_r-5,1);
                        Zpatch(1:5,i_patch)   = [z_pre(jj) z_here(jj) z_here(jj+1) z_pre(jj+1) z_pre(jj)]';
                        Zpatch(6:end,i_patch) = z_pre(jj)*ones(n_r-5,1);
                        i_patch = i_patch+1;

                    end

                    x_pre    = x_here;
                    y_pre    = y_here;
                    z_pre    = z_here;

                end

                Xpatch(:,i_patch) = x_here';
                Ypatch(:,i_patch) = y_here';
                Zpatch(:,i_patch) = z_here';

                patch(Xpatch,Ypatch,Zpatch,color,'EdgeColor','none','FaceAlpha',1);
            else
                CustomShapePlot(g_here);
            end
            gf     = Tr.VLinks(Tr.LinkIndex(i)).gf;
            g_here = g_here*gf;
            Ad_gf_inv = dinamico_Adjoint(ginv(gf));
            J_here    = Ad_gf_inv*J_here;
            
        end
        end
        
        dof_start = dof_start+dof_here;
        
            %=============================================================================
        for j=1:(Tr.VLinks(Tr.LinkIndex(i)).npie)-1 % for each piece
            
            dof_here   = Tr.CVTwists{i}(j+1).dof;
            Type       = Tr.CVTwists{i}(j+1).Type;
            q_here     = q(dof_start:dof_start+dof_here-1);
            xi_starfn  = Tr.CVTwists{i}(j+1).xi_starfn;
            gi         = Tr.VLinks(Tr.LinkIndex(i)).gi{j};
            Bdof      = Tr.CVTwists{i}(j+1).Bdof;
            Bodr      = Tr.CVTwists{i}(j+1).Bodr;
            Bh         = Tr.CVTwists{i}(j+1).Bh;
            ld        = Tr.VLinks(Tr.LinkIndex(i)).lp{j};
            g_here     = g_here*gi;
            Ad_gi_inv = dinamico_Adjoint(ginv(gi));
            J_here    = Ad_gi_inv*J_here;
               
            Xs          = linspace(0,1,n_l);
            color       = Tr.VLinks(Tr.LinkIndex(i)).color;
            H           = Xs(2)-Xs(1);
            Z           = (1/2)*H;          % Zanna quadrature coefficient

            Xpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
            Ypatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
            Zpatch  = zeros(n_r,(n_r-1)*(n_l-2)+2);
            i_patch = 1;
            
            if Tr.VLinks(Tr.LinkIndex(i)).CS=='C'
                
                r_fn  = Tr.VLinks(Tr.LinkIndex(i)).r{j};
                r     = r_fn(0);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = r*sin(theta);
                z     = r*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];
            elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='R'
                h_fn = Tr.VLinks(Tr.LinkIndex(i)).h{j};
                w_fn = Tr.VLinks(Tr.LinkIndex(i)).w{j};
                h    = h_fn(0);
                w    = w_fn(0);
                x    = [0 0 0 0 0];
                y    = [h/2 -h/2 -h/2 h/2 h/2];
                z    = [w/2 w/2 -w/2 -w/2 w/2];
                pos  = [x;y;z;ones(1,5)];
                
            elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='E'
                a_fn  = Tr.VLinks(Tr.LinkIndex(i)).a{j};
                b_fn  = Tr.VLinks(Tr.LinkIndex(i)).b{j};
                a     = a_fn(0);
                b     = b_fn(0);
                theta = linspace(0,2*pi,n_r);
                x     = zeros(1,n_r);
                y     = a*sin(theta);
                z     = b*cos(theta);
                pos   = [x;y;z;ones(1,n_r)];
            end

            pos_here = g_here*pos;
            x_here   = pos_here(1,:);
            y_here   = pos_here(2,:);
            z_here   = pos_here(3,:);

            Xpatch(:,i_patch) = x_here';
            Ypatch(:,i_patch) = y_here';
            Zpatch(:,i_patch) = z_here';
            i_patch           = i_patch+1;

            x_pre = x_here;
            y_pre = y_here;
            z_pre = z_here;
            
            Lscale = ld;
            
            for ii=1:n_l-1
                
                if Tr.VLinks(Tr.LinkIndex(i)).CS=='C'
                    r     = r_fn(Xs(ii+1));
                    theta = linspace(0,2*pi,n_r);
                    x     = zeros(1,n_r);
                    y     = r*sin(theta);
                    z     = r*cos(theta);
                    pos   = [x;y;z;ones(1,n_r)];
                elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='R'

                    h   = h_fn(Xs(ii+1));
                    w   = w_fn(Xs(ii+1));
                    x   = [0 0 0 0 0];
                    y   = [h/2 -h/2 -h/2 h/2 h/2];
                    z   = [w/2 w/2 -w/2 -w/2 w/2];
                    pos = [x;y;z;ones(1,5)];
                elseif Tr.VLinks(Tr.LinkIndex(i)).CS=='E'

                    a     = a_fn(Xs(ii+1));
                    b     = b_fn(Xs(ii+1));
                    theta = linspace(0,2*pi,n_r);
                    x     = zeros(1,n_r);
                    y     = a*sin(theta);
                    z     = b*cos(theta);
                    pos   = [x;y;z;ones(1,n_r)];
                end
                
                X   = Xs(ii);
                X_Z = X+Z;
                
                xi_Zhere  = xi_starfn(X_Z);
                
                if ~isempty(q_here)
                    if strcmp(Type,'FEM Like')
                        SubClass = Tr.CVTwists{i}(j+1).SubClass;
                        B_Zhere  = diag([1/Lscale 1/Lscale 1/Lscale 1 1 1])*Bh(X_Z,Bdof,Bodr,SubClass);
                    elseif strcmp(Type,'Custom Independent')
                        B_Zhere  = diag([1/Lscale 1/Lscale 1/Lscale 1 1 1])*Bh(X_Z);
                    else
                        B_Zhere  = diag([1/Lscale 1/Lscale 1/Lscale 1 1 1])*Bh(X_Z,Bdof,Bodr);
                    end
                end

                xi_Zhere  = xi_Zhere+B_Zhere*q_here;
                Gamma_here  = H*Lscale*xi_Zhere;

                [gh,Tg]    = variable_expmap_gTg_mex(Gamma_here);
                TBGamma_here                                    = zeros(6,Tr.ndof);
                TBGamma_here(:,dof_start:dof_start+dof_here-1)  = Tg*H*Lscale*B_Zhere;

                g_here     = g_here*gh;
                Ad_gh_inv  = dinamico_Adjoint(ginv(gh));
                J_here     = Ad_gh_inv*(J_here+TBGamma_here); %full
                
                pos_here = g_here*pos;
                x_here   = pos_here(1,:);
                y_here   = pos_here(2,:);
                z_here   = pos_here(3,:);


                for jj=1:n_r-1

                    Xpatch(1:5,i_patch)   = [x_pre(jj) x_here(jj) x_here(jj+1) x_pre(jj+1) x_pre(jj)]';
                    Xpatch(6:end,i_patch) = x_pre(jj)*ones(n_r-5,1);
                    Ypatch(1:5,i_patch)   = [y_pre(jj) y_here(jj) y_here(jj+1) y_pre(jj+1) y_pre(jj)]';
                    Ypatch(6:end,i_patch) = y_pre(jj)*ones(n_r-5,1);
                    Zpatch(1:5,i_patch)   = [z_pre(jj) z_here(jj) z_here(jj+1) z_pre(jj+1) z_pre(jj)]';
                    Zpatch(6:end,i_patch) = z_pre(jj)*ones(n_r-5,1);
                    i_patch = i_patch+1;

                end
                
                x_pre = x_here;
                y_pre = y_here;
                z_pre = z_here;
                
            end

            Xpatch(:,i_patch) = x_here';
            Ypatch(:,i_patch) = y_here';
            Zpatch(:,i_patch) = z_here';

            patch(Xpatch,Ypatch,Zpatch,color,'EdgeColor','none');
            
%             g_here(1:3,4)  = g_here(1:3,4)*Lscale;
            %updating g, Jacobian, Jacobian_dot and eta at X=L
            gf     = Tr.VLinks(Tr.LinkIndex(i)).gf{j};
            g_here = g_here*gf;
            Ad_gf_inv = dinamico_Adjoint(ginv(gf));
            J_here    = Ad_gf_inv*J_here;
            
            dof_start = dof_start+dof_here;
            
        end
        g_Ltip((i-1)*4+1:i*4,:) = g_here;
        J_Ltip((i-1)*6+1:i*6,:) = J_here;

    end
    ptip(:,it) = g_here(1:3,4);
    gtip(4*(it-1)+1:4*it,:) = g_here;
    eta_tip = J_here*qd;

    p1     = g_here*[0 0 0 1]';
    vdir   = eta_tip(4:6)/norm(eta_tip(4:6));
    vdir   = g_here(1:3,1:3)*vdir; %to global
    % vdir   = vdir(1:3);
    
    scale = 0.03;
    p1 = p1(1:3);
    if norm(eta_tip(4:6))>0
        p2 = p1+vdir*scale*norm(eta_tip(4:6));
        hv1 = arrow3(p1',p2','r2',0.5);
    end
    % i1=1;
    if it<FrameRate
        i1=1;
    else
        i1 = it-FrameRate+1;
    end
    h2 = plot3(ptip(1,i1:it),ptip(2,i1:it),ptip(3,i1:it),'k','LineWidth',1.5);
    it = it+1;
    
    scatter3(xyz_gt_pos(1), xyz_gt_pos(2), xyz_gt_pos(3), 100, 'r', 'filled'); % 红色标记目标点
    scatter3(xyz_pred_pos(1), xyz_pred_pos(2), xyz_pred_pos(3), 100, 'b', 'filled'); % 红色标记目标点
    frame = getframe(gcf);
    writeVideo(v,frame);
end

close(v);
end
% 
% answer = questdlg('Play output video in MATLAB?','Grapical Output', ...
% 	'Yes','No','Yes');
% 
% if strcmp('Yes',answer)
%     implay('.\Dynamics.avi')
% end
