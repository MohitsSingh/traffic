%% BATCHTRAFFIC
% Run the traffic simulation in a loop
% Version: 2015jun11

%% Parameters and settings

tic
disp('Setting parameters...')
% General code stuff
doplot = 0; % Toggle whether to plot
nsims = 10; % Number of simulations
nonparametric = true;

% Boring initialization stuff
maxtime = 10; % Number of seconds to simulate
dx = 0.1; % Spatial step for calculating look-up table
dt = 0.02; % Timestep size in seconds
npts = maxtime/dt;
time = dt:dt:maxtime;

% Model parameters
whichmodel = 5;
useawake = 1;
tracklength = 300; % Track length in meters
tau = round(1.2/dt); % Reaction time (in number of timesteps)
randvel = 0; % Randomness in initial velocity in m/s
randpos = 0; % Randomness in initial positions in m

% User-set variables
ndriverslist = 20; % Number of drivers
defaultvlist = [5, 10, 20, 50, 100]; % Default velocities in m/s
ninputs = 0;
inputs = struct();
for d = 1:length(ndriverslist)
    for v = 1:length(defaultvlist)
        ninputs = ninputs+1;
        inputs(ninputs).ndrivers = ndriverslist(d);
        inputs(ninputs).defaultv = defaultvlist(v);
    end
end

% inputs.fracinattentive = 0.5; % Fraction of drivers who are inattentive




%% Markov process parameters
if useawake == 1
    tausleep = 6; % Time constant for losing alertness
    tauwake =  3; % Time constant for regaining alertness
else
    tausleep = 1e9; % Never sleep
    tauwake = 1e-9; % Always wake, watching
end


%% Define outputs
nsimzeros = zeros(ninputs, nsims);
outputs = struct();
outputs.collisions = nsimzeros;
outputs.meanvelocity = nsimzeros;
outputs.headwaycv = nsimzeros;


%% Define traffic model
switch whichmodel
    
    %% Bando 98 Model
    case 1
        % Set parameters
        C(1) = 16.8;
        C(2) = 0.086;
        C(3) = C(2)*-25;
        C(4) = C(1)*0.913;
        vmax = C(1)+C(4);
        alpha = 0.5; % Responsivity of drivers
        % Define optimal velocity function
        dvdt = @(perchead,percv,percvdiff) alpha*dt*(max(0,C(1)*tanh(C(2)*perchead+C(3))+C(4))-percv); % Taken from paper

    
    %% Orosz Model
    case 2
        disp('CHECK')
        % Set parameters
        hstop = 7;
        vmax = 32;
        alpha = 0.5; % Responsivity of drivers
        % Define optimal velocity function
        dvdt = @(perchead,percv,percvdiff) alpha*dt*(max(0,vmax*(perchead/hstop - 1).^3./(1 + (perchead/hstop - 1).^3))-percv);
        
    %% Human Driver Model
    case 3
        w = 0; % Initialization of the Wiener process
        corrt = 20; % Correlation time in s
        na = 1; % Number of anticipated vehicles (look-ahead)
        relerror = 0.05; % Relative error in estimating headway
        rc = 0.01; % Inverse "TTC" error in /s
        vmax = 128/3.6; % Maximum velocity in m/s
        T = 1.1; % ?? in s
        a = 1; % ??? in m/s^2
        b = 1.5; % ??? m/s^2
        s0 = 2; % Minimum length in m
        wiener = @(w) exp(-dt/corrt)*w+sqrt(2*dt/corrt)*randn(size(w));
        
        
    %% Generalized Force Model
    case 4
        alpha = 1/2.45; % Responsivity rate in /s
        vmax = 16.98; % ?? Maximum velocity ish?
        d = 1.38; % Minimum vehicle distance in m
        T = 0.74; % Safe time headway in s
        alphaprime = 1/0.77; % Braking rate in s
        R = 5.59; % ?? in m
        Rprime = 98.78; % Braking interaction distance in m
        
        s = @(percv) d+T*percv;
        V = @(perchead,percv) vmax*(1-exp(-(perchead-s(percv))/R));
        dvdt = @(perchead,percv,percvdiff) alpha*dt*(vmax - percv + V(perchead,percv)-vmax) - alphaprime*dt*max(0,percvdiff)*exp(-(perchead-s(percv))/Rprime);
        
        
    %% Underwood Model
    case 5
        vmax = 24.40; % Maximum velocity in m/s
        hm = 12.92; % Headway inflection point in m
        alpha = 0.88; % Responsivity in /s
        dvdt = @(perchead,percv,percvdiff) alpha*dt*(vmax*exp(-2*hm./perchead)-percv);
    
    otherwise
        error('OMG 5 options and you still can''t decide. laaame.')
end




%% Loop over inputs and simulations
disp('Running simulations...')
loopcount = 0;
for i = 1:ninputs % Loop over number of drivers, for example...not implemented, for now
    ndrivers = inputs(i).ndrivers;
    defaultv = inputs(i).defaultv;
    for s = 1:nsims
        loopcount = loopcount+1;
        fprintf('  %i of %i\n', loopcount, ninputs*nsims)
        
        %% Set up initial conditions
        
        drivers = struct;
        for d = 1:ndrivers
            drivers(d).x = (d/ndrivers)*tracklength+randpos*randn; % Initial position
            drivers(d).v = defaultv+randvel*randn; % Initial velocity
            drivers(d).alert = 1; % Alert vs. fucked
            drivers(d).realhead = tracklength/ndrivers; % Amount of distance to car in front
            drivers(d).perchead = tracklength/ndrivers; % Amount of distance to car in front
            drivers(d).percv = 0; % Perceived velocity
            drivers(d).percvdiff = 0; % Perceived velocity difference
            drivers(d).infront = mod(d,ndrivers)+1; % Index of driver in front
            drivers(d).behind = mod(d-2,ndrivers)+1; % Index of driver behind
        end
        
        
        %% Run simulation
        
        % Initialize data matrices
        ncollisions = 0;
        matrix = zeros(ndrivers, npts);
        positions = matrix;
        velocities = matrix;
        realheads = matrix;
        percheads = matrix;
        percvs = matrix;
        percvdiffs = matrix;
        alert = matrix;
        collisions = matrix;
        
        % Solve difference equations across all time points
        for t = 1:npts
            
            % Solve motion of each driver at timepoint 'd'
            for d = 1:ndrivers
                if t>tau+1 % Update velocity
                    sr = rand(1); % Random variable for determing transitions of Markov process
                    if drivers(d).alert == 1
                        if sr < dt/tausleep
                            drivers(d).alert = 0;
                        end
                    else
                        if sr < dt/tauwake
                            drivers(d).alert = 1;
                        end
                    end
                    
                    if drivers(d).alert == 1,
                        drivers(d).percv = velocities(d,t-tau-1); % Update perceived velocity if driver is alert
                        drivers(d).perchead = mod(diff(positions([d drivers(d).infront],t-tau-1)),tracklength); % Update perceived headway if driver is alert
                        drivers(d).percvdiff = diff(velocities([drivers(d).infront d],t-tau-1)); % Update perceived velocity difference if driver is alert
                    end
                    drivers(d).v = drivers(d).v+dvdt(drivers(d).perchead, drivers(d).percv, drivers(d).percvdiff); % Update velocity as a function of perceived headway
                end
                drivers(d).x = mod(drivers(d).x+drivers(d).v*dt,tracklength); % Update position
                
                positions(d,t) = drivers(d).x; % Save current position
                velocities(d,t) = drivers(d).v; % Save current position
                percvs(d,t) = drivers(d).percv; % Save current position
                realheads(d,t) = drivers(d).realhead; % Save current actual headway
                percheads(d,t) = drivers(d).perchead; % Save current perceived head
                percvdiffs(d,t) = drivers(d).v; % Save current position
                alert(d,t) = drivers(d).alert; % Save state
            end
            
            % Use actual headway to check for drivers overtaking one another; treat these events as collisions
            for d = 1:ndrivers
                % Update positions
                if t>1
                    
                    drivers(d).realhead = mod(diff(positions([d drivers(d).infront],t-1)),tracklength);
                    % Check if collision occurred. If it did, reorder drivers and
                    % alter velocities of collided vehicles
                    count = 0;
                    while drivers(d).realhead>tracklength-2*dt*vmax
                        count = count+1;
                        ncollisions = ncollisions+1;
                        me = d;
                        behind = drivers(me).behind;
                        infront1 = drivers(me).infront;
                        infront2 = drivers(infront1).infront;
                        drivers(behind).infront = infront1;
                        drivers(me).infront = infront2;
                        drivers(me).behind = infront1;
                        drivers(infront1).infront = me;
                        drivers(infront1).behind = behind;
                        drivers(infront2).behind = me;
                        drivers(d).realhead = mod(diff(positions([d drivers(d).infront],t-1)),tracklength);
                        drivers(d).v = 0;
                        drivers(drivers(d).behind).v = 0;
                        collisions(d,t) = 1;
                    end
                end
            end
            
        end
        
        outputs.collisions(i,s) = ncollisions;
        outputs.meanvelocity(i,s) = mean(velocities(:));
        outputs.headwaycv(i,s) = std(percheads(:))/mean(percheads(:));
        
    end
end

%% Data analysis
outputnames = {'Number of collisions', 'Mean velocity', 'Headway CV'};
stats = struct();
outputfields = fieldnames(outputs);
noutputs = length(outputfields);
for i = 1:noutputs
    f = outputfields{i};
    stats.(f) = struct();
    if nonparametric
        quantiles = quantile(outputs.(f), [0.25, 0.5, 0.75], 2);
        stats.(f).best = quantiles(:, 2);
        stats.(f).low = quantiles(:, 1);
        stats.(f).high = quantiles(:, 3);
    else
        stats.(f).best = mean(outputs.(f), 2);
        stats.(f).std = std(outputs.(f), 2);
        stats.(f).low = stats.(f).best - stats.(f).std;
        stats.(f).high = stats.(f).best + stats.(f).std;
    end
end


%% Plotting
figure()
for i=1:noutputs
    subplot(2,2,i); hold on
    f = outputfields{i};
    errorbar(defaultvlist, stats.(f).best, stats.(f).best-stats.(f).low, stats.(f).high-stats.(f).best)
    scatter(defaultvlist, stats.(f).best)
    title(outputnames{i})
end




disp('Done.')