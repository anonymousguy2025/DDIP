% Initial conditions
clc
load('SerialRobot','S1','q0','u0')
S1.PlotParameters.ClosePrevious=false;
qu0 = q0;
qu0(S1.i_jactq) = u0;

n = 1000;
uq = pi/4-pi/2*rand(S1.nact,n);

%% method 1 SoRoSim Algorithm
options = optimoptions('fsolve','Display','off','Algorithm','trust-region-dogleg');%set display to 'final' to see if solved. 'off' for none. 'trust-region-dogleg' (default), 'trust-region', and 'levenberg-marquardt'

qu_soro = zeros(S1.ndof,n);
q_soro = zeros(S1.ndof,n);

tic
for i=1:n
qu_soro(:,i) = fsolve(@(qu) Equilibrium_SoRoSim(S1,qu,uq(:,i)),qu0,options); %Residue for only residue, ResidueJacobian for Residue and Jacobian
q_soro(:,i) = qu_soro(:,i);
q_soro(S1.i_jactq,i) = uq(:,i);

end
toc
save('StaticsSolution_SoRoSim.mat','uq','qu_soro','q_soro')

figure(1)

ipick = randi([1, 1000], 1, 5);
for i=ipick
    if i==ipick(1)
        S1.PlotParameters.Light=true;
        S1.plotq(q_soro(:,i))
        S1.PlotParameters.Light=false;
    else
        S1.plotq(q_soro(:,i))
    end
end
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



%% method 3, RNEA using Jacobian of Residue
options = optimoptions('fsolve','Jacobian','on','Display','off','Algorithm','trust-region-dogleg');%set display to 'final' to see if solved. 'off' for none. 'trust-region-dogleg' (default), 'trust-region', and 'levenberg-marquardt'


qu_J = zeros(S1.ndof,n);
q_J = zeros(S1.ndof,n);

tic
for i=1:n

qu_J(:,i) = fsolve(@(q) StaticsResidueJacobian(S1,q,uq(:,i)),q0,options); %Residue for only residue, ResidueJacobian for Residue and Jacobian
q_J(:,i) = qu_J(:,i);
q_J(S1.i_jactq,i) = uq(:,i);

end
toc

save('StaticsSolution_SoRoSimwithJacobian.mat','uq','qu_J','q_J')

figure(2)

for i=ipick
    if i==ipick(1)
        S1.PlotParameters.Light=true;
        S1.plotq(q_J(:,i))
        S1.PlotParameters.Light=false;
    else
        S1.plotq(q_J(:,i))
    end
end
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