
%{
Types of Neuron Models used
1 - regular spiking (excitatory)
2 - fast spiking (inhibitory)
3 - low-threshold spiking (inhibitory)
4 - chattering neurons (excitatory)
%}

warning('off')
% Number of Neurons
prompt = "Number of Neurons to Simulate: ";
x = input(prompt);
num_neu = x;

% Initialization of Neuron types
neu_list = zeros(num_neu,1);
M = size(neu_list,1);
idx_neu_list = randperm(numel(neu_list),round(numel(neu_list))); 

prompt = "Percentage of regular spiking neurons desired (input a number between 1-100): ";
x1 = input(prompt);
x1 = round(x1*num_neu/100);
neu_list(idx_neu_list(1:x1)) = 1;

prompt = "Percentage of fast spiking desired (input a number between 1-100): ";
x2 = input(prompt);
x2 = round(x2*num_neu/100);
neu_list(idx_neu_list(x1+1:x1+x2)) = 2;

prompt = "Percentage of low-threshold spiking neurons desired (input a number between 1-100). Percentage of chattering neurons will be set automatically after this: ";
x3 = input(prompt);
x3 = round(x3*num_neu/100);
neu_list(idx_neu_list(x1+x2+1:x1+x2+x3)) = 3;

neu_list(idx_neu_list(x1+x2+x3+1:end)) = 4;
neu_list


% Simulation parameters
prompt = "Set time step in seconds: ";
dt = input(prompt);
prompt = "Set total simulation time in seconds: ";
t_tot = input(prompt);
T = t_tot/dt;

% Different current initializations for applied current, synapse noise,
% synpatical connections, and total contributions of all the above
Iapp = zeros(T,M);
Isyn_noise = zeros(T,M);
Isyn = zeros(T,M);
Itot = zeros(T,M);

% Initializing Membrane Potential-related values
v = zeros(T,M); % Membrane Potential
v(1,:) = -70; % Initial Membrane Potential
u = zeros(T,M); % Membrane recovery potential
u(1,:) = -14; % Initial Membrane recovery potential

% Synapse model paramters
tg = 10; % time constant for first-order linear kinetics
frate = 0.002; % Rate of Possion Process

prompt = "How many synapses do you want to define for noise modelling? ";
nsyn = input(prompt);
% nsyn = 100; % Number of input synapses for purpose of noise modelling 
% (can be tuned)
r_nn = rand(T,nsyn,M); % Random numbers generated for Poisson Process
p_nn = zeros(nsyn,M); % Variable to update conductances based on Poisson Process
gin_nn = zeros(nsyn,M); % Conductances for Poisson Process
E_nn = zeros(nsyn,M); % Used as a variable in the Possion Process modelling
% of synapses
win_nn = zeros(nsyn,M); % Weight's for Possion Process
% Select 15% random indices
idx_nn = randperm(numel(win_nn),round(numel(win_nn)*15/100)); 
% 0.1 weight assigned to 15% of all noisy connections
win_nn(idx_nn) = 0.1;


% Recurrent Parameters
gin_ij = zeros(M,M); % Conductances for neuron i to neuron j connections
E_ij = zeros(M,M); % Variable to help model synaptic connections
win_ij = zeros(M,M); % Weights for connections from neuron i to neuron j
% Select 15% random indices
idx_ij = randperm(numel(win_ij),round(numel(win_ij)*15/100));
% 15% of all recurrent connections have weights pulled out randomly from a
% Gaussian Distribution
win_ij(idx_ij) = gamrnd(4,0.007,1,round(numel(win_ij)*15/100));
% Make the matrix sparse, thereby allowing faster computation
win_ij = sparse(win_ij);

prompt = "How many neuron connections do you want to define? ";
nc_num = input(prompt);

for n = 1:nc_num
    prompt = "Neurons to connect to each other. Enter Neuron 1 index: ";
    win1 = input(prompt);
    prompt = "Neurons to connect to each other. Enter Neuron 2 index: ";
    win2 = input(prompt);    
    win_ij(win2, win1) = 1;
end

for t = 1:T-1
    for ii = 1:M
        ntype = neu_list(ii);
        % regular spiking
        if ntype == 1
            a = 0.02; % time scale of u
            b = 0.2; % sensitivity of recovery variable
            % after-spike reset value of v
            c = -65; % Deep voltage reset
            % after-spike reset of recovery variable
            d = 8; % large after-spike jump of u
            vthresh = 35; % threshold voltage
        end
        % fast spiking
        if ntype == 2
            a = 0.1; % fast recovery
            b = 0.2;
            c = -65;
            d = 2;
            vthresh = 35;
            % If inhibitory neuron, set Ej = -85
            E_ij(:,ii) = - 85;
        end
        % low-threshold spiking
        if ntype == 3
            a = 0.02; % fast recovery
            b = 0.25; % low firing thresholds
            c = -65;
            d = 2;
            vthresh = 35;
            % If inhibitory neuron, set Ej = -85
            E_ij(:,ii) = - 85;
        end
        % chattering neurons
        if ntype == 4
            a = 0.02;
            b = 0.2;
            c = -50; % very high voltage reset  
            d = 2; % moderate after-spike jump of u
            vthresh = 35;
        end
        for s = 1:nsyn
            if (t*dt > 150) && (t*dt < 800)
                Iapp(t,ii) = 7;
                if (r_nn(t,s,ii) <= frate*dt)
                    p_nn(s,ii) = 1;
                else
                    p_nn(s,ii) = 0;
                end
            else
                Iapp(t,ii) = 0;
                p_nn(s,ii) = 0;
            end
        end
        gin_nn(:,ii) = gin_nn(:,ii) + p_nn(:,ii);
        Isyn_noise(t,ii) = win_nn(:,ii).'*(gin_nn(:,ii).*E_nn(:,ii)) - ((win_nn(:,ii).')*gin_nn(:,ii))*v(t,ii);
        gin_nn = (1-dt/tg)*gin_nn;
        if v(t,ii) >= vthresh
            gin_ij(:,ii) = gin_ij(:,ii) + 1;
        end
%         if mod(ii,2) == 0
%             win_ij(ii-1,ii) = 1;
%             win_ij(ii,ii-1) = 1;
%         end
        Isyn(t,ii) = win_ij(:,ii).'*(gin_ij(:,ii).*E_ij(:,ii)) - ((win_ij(:,ii).')*gin_ij(:,ii))*v(t,ii);
        gin_ij = (1-dt/tg)*gin_ij;
        Itot(t,ii) = Iapp(t,ii) + Isyn_noise(t,ii) + Isyn(t,ii);
        if v(t,ii) < vthresh
            v(t+1,ii) = v(t,ii) + dt*(0.04*v(t,ii)^2 + 5*v(t,ii) + 140 - u(t,ii) + Itot(t,ii));
            u(t+1,ii) = u(t,ii) + dt*a*(b*v(t,ii)-u(t,ii));
        else
            v(t+1,ii) = c;
            u(t+1,ii) = u(t,ii) + d;
            v(t,ii) = 35;
            % Plotting spike raster
            if (ntype == 1 || ntype == 4)
                scatter(t,ii,15,"filled","black")
                hold on
                xlabel('Time');
                ylabel('Neuron Idx');
            else
                scatter(t,ii,15,"filled","red")
                hold on
                xlabel('Time');
                ylabel('Neuron Idx');
            end
        end
    end
end

figure
% f = figure('WindowState','maximized');
% pause(1);
tt = tiledlayout('flow');
tt.TileSpacing = 'compact';
tt.Padding = 'compact';

for ii = 1:M
    nexttile()
    plot(1:T, v(:,ii));
    grid on
    hold on
    plot(1:T, Iapp(:,ii));
    xlabel('Time');
    ylabel(['v N',num2str(ii)]);
end

figure
% f = figure('WindowState','maximized');
% pause(1);
for ii = 1:M
    nexttile()
    plot(1:T, Itot(:,ii));
    grid on
    hold on
    xlabel('Time');
    ylabel(['Itot N', num2str(ii)]);
end