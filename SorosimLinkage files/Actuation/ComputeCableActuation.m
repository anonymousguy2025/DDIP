%Function to compute cable actuation when cable is fully inside (14.06.2021)
function Bq = ComputeCableActuation(Tr,dc,dcp,Sdiv,Ediv,q)

Bq        = zeros(Tr.ndof,1);
dof_start = 1;

for i=1:Tr.N

    dof_start = dof_start+Tr.CVTwists{i}(1).dof;
    
    for j=1:Tr.VLinks(Tr.LinkIndex(i)).npie-1
        
        VTwists  = Tr.CVTwists{i};
        dof_here = VTwists(j+1).dof;

        if j>=Sdiv(i)&&j<=Ediv(i)

            dcd  = dc{i}{j};
            dcpd = dcp{i}{j};
            
            q_here  = q(dof_start:dof_start+dof_here-1);
            xi_star = VTwists(j+1).xi_star;
            Ws      = VTwists(j+1).Ws;
            nip     = VTwists(j+1).nip;

            for k=1:nip
                if Ws(k)>0
                    dc_here     = dcd(:,k);
                    dcp_here    = dcpd(:,k);
                    xi_here = xi_star(6*(k-1)+1:6*k,1);

                    B_here = VTwists(j+1).B((k-1)*6+1:k*6,:);
                    if dof_here>0
                        xi_here    = B_here*q_here+xi_here;
                    end

                    xihat_here123  = [0 -xi_here(3) xi_here(2) xi_here(4);xi_here(3) 0 -xi_here(1) xi_here(5);-xi_here(2) xi_here(1) 0 xi_here(6)];%4th row is avoided to speedup calculation
                    tang        = xihat_here123*[dc_here;1]+dcp_here;
                    tang        = tang/norm(tang);
                    Btau        = [[0 -dc_here(3) dc_here(2);dc_here(3) 0 -dc_here(1);-dc_here(2) dc_here(1) 0]*tang;tang];

                    Bq(dof_start:dof_start+dof_here-1)  = Bq(dof_start:dof_start+dof_here-1)+Ws(k)*B_here'*Btau; %scaled back for addition
                end
            end
        end
        dof_start = dof_start+dof_here;
    end 
end

end

