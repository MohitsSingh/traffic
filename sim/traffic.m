%% TRAFFIC
% Our totally awesome traffic model.
% Version: 2012nov04

%% Parameters and settings

tic
disp('Setting parameters...')
% General code stuff
rng(239487) % Random number seed
doplot=1; % Toggle whether to plot

% Boring initialization stuff
ncollisions=0; % Initialize number of collisions
maxtime=50; % Number of seconds to simulate
dx=0.1; % Spatial step for calculating look-up table
dt=0.01; % Timestep size in seconds
npts=maxtime/dt;
time=dt:dt:maxtime;

% Model parameters
whichmodel=5;
useawake=1;
ndrivers=100; % Number of drivers
tracklength=3000; % Track length in meters
defaultv=20; % Default velocity in m/s
tau=round(0.0/dt); % Reaction time (in number of timesteps)
randvel=5; % Randomness in initial velocity in m/s
randpos=5; % Randomness in initial positions in m


%% Markov process parameters
if useawake==1
    tausleep = 6; % Time constant for losing alertness
    tauwake =  3; % Time constant for regaining alertness
else
    tausleep = 1e9;
    tauwake = 1e-9;
end


%% Set up initial conditions

drivers=struct;
for d=1:ndrivers
    drivers(d).x=(d/ndrivers)*tracklength+randpos*randn; % Initial position
    drivers(d).v=defaultv+randvel*randn; % Initial velocity
    drivers(d).alert=1; % Alert vs. fucked
    drivers(d).realhead=tracklength/ndrivers; % Amount of distance to car in front
    drivers(d).perchead=tracklength/ndrivers; % Amount of distance to car in front
    drivers(d).percv=0; % Perceived velocity
    drivers(d).percvdiff=0; % Perceived velocity difference
    drivers(d).infront=mod(d,ndrivers)+1; % Index of driver in front
    drivers(d).behind=mod(d-2,ndrivers)+1; % Index of driver behind
end



%% Define traffic model
switch whichmodel
    
    %% Bando 98 Model
    case 1
        % Set parameters
        C(1)=16.8;
        C(2)=0.086;
        C(3)=C(2)*-25;
        C(4)=C(1)*0.913;
        vmax=C(1)+C(4);
        alpha=0.5; % Responsivity of drivers
        % Define optimal velocity function
        dvdt=@(perchead,percv,percvdiff) alpha*dt*(max(0,C(1)*tanh(C(2)*perchead+C(3))+C(4))-percv); % Taken from paper

    
    %% Orosz Model
    case 2
        disp('CHECK')
        % Set parameters
        hstop=7;
        vmax=32;
        alpha=0.5; % Responsivity of drivers
        % Define optimal velocity function
        dvdt=@(perchead,percv,percvdiff) alpha*dt*(max(0,vmax*(perchead/hstop - 1).^3./(1 + (perchead/hstop - 1).^3))-percv);
        
    %% Human Driver Model
    case 3
        w=0; % Initialization of the Wiener process
        corrt=20; % Correlation time in s
        na=1; % Number of anticipated vehicles (look-ahead)
        relerror=0.05; % Relative error in estimating headway
        rc=0.01; % Inverse "TTC" error in /s
        vmax=128/3.6; % Maximum velocity in m/s
        T=1.1; % ?? in s
        a=1; % ??? in m/s^2
        b=1.5; % ??? m/s^2
        s0=2; % Minimum length in m
        wiener=@(w) exp(-dt/corrt)*w+sqrt(2*dt/corrt)*randn(size(w));
        
        
    %% Generalized Force Model
    case 4
        alpha=1/2.45; % Responsivity rate in /s
        vmax=16.98; % ?? Maximum velocity ish?
        d=1.38; % Minimum vehicle distance in m
        T=0.74; % Safe time headway in s
        alphaprime=1/0.77; % Braking rate in s
        R=5.59; % ?? in m
        Rprime=98.78; % Braking interaction distance in m
        
        s=@(percv) d+T*percv;
        V=@(perchead,percv) vmax*(1-exp(-(perchead-s(percv))/R));
        dvdt=@(perchead,percv,percvdiff) alpha*dt*(vmax - percv + V(perchead,percv)-vmax) - alphaprime*dt*max(0,percvdiff)*exp(-(perchead-s(percv))/Rprime);
        
        
    %% Underwood Model
    case 5
        vmax=24.40; % Maximum velocity in m/s
        hm=12.92; % Headway inflection point in m
        alpha=0.88; % Responsivity in /s
        dvdt=@(perchead,percv,percvdiff) alpha*dt*(vmax*exp(-2*hm./perchead)-percv);
    
    otherwise
        error('OMG 5 options and you still can''t decide. laaame.')
end



%% Run simulation

% Initialize data matrices
disp('Running simulation...')
matrix=zeros(ndrivers,npts);
positions=matrix;
velocities=matrix;
realheads=matrix;
percheads=matrix;
percvs=matrix;
percvdiffs=matrix;
alert=matrix;
collisions=matrix;

% Solve difference equations across all time points
for t=1:npts
    
    % Solve motion of each driver at timepoint 'd'
    for d=1:ndrivers
        if t>tau+1 % Update velocity
            sr = rand(1); % Random variable for determing transitions of Markov process
            if drivers(d).alert==1
                if sr < dt/tausleep
                    drivers(d).alert=0;
                end
            else
                if sr < dt/tauwake
                    drivers(d).alert=1;
                end
            end
            
            if drivers(d).alert == 1,
                drivers(d).percv=velocities(d,t-tau-1); % Update perceived velocity if driver is alert
                drivers(d).perchead=mod(diff(positions([d drivers(d).infront],t-tau-1)),tracklength); % Update perceived headway if driver is alert
                drivers(d).percvdiff=diff(velocities([drivers(d).infront d],t-tau-1)); % Update perceived velocity difference if driver is alert
            end
            drivers(d).v=drivers(d).v+dvdt(drivers(d).perchead, drivers(d).percv, drivers(d).percvdiff); % Update velocity as a function of perceived headway
        end
        drivers(d).x=mod(drivers(d).x+drivers(d).v*dt,tracklength); % Update position
        
        positions(d,t)=drivers(d).x; % Save current position
        velocities(d,t)=drivers(d).v; % Save current position
        percvs(d,t)=drivers(d).percv; % Save current position
        realheads(d,t)=drivers(d).realhead; % Save current actual headway
        percheads(d,t)=drivers(d).perchead; % Save current perceived head
        percvdiffs(d,t)=drivers(d).v; % Save current position
        alert(d,t)=drivers(d).alert; % Save state
    end
    
    % Use actual headway to check for drivers overtaking one another; treat these events as collisions
    for d=1:ndrivers
        % Update positions
        if t>1
  
            drivers(d).realhead=mod(diff(positions([d drivers(d).infront],t-1)),tracklength);
            % Check if collision occurred. If it did, reorder drivers and
            % alter velocities of collided vehicles
            count=0;
            while drivers(d).realhead>tracklength-2*dt*vmax
                count=count+1;
                ncollisions=ncollisions+1;
                me=d;
                behind=drivers(me).behind;
                infront1=drivers(me).infront;
                infront2=drivers(infront1).infront;
                drivers(behind).infront=infront1;
                drivers(me).infront=infront2;
                drivers(me).behind=infront1;
                drivers(infront1).infront=me;
                drivers(infront1).behind=behind;
                drivers(infront2).behind=me;
                drivers(d).realhead=mod(diff(positions([d drivers(d).infront],t-1)),tracklength);
                drivers(d).v=0;
                drivers(drivers(d).behind).v=0;
                collisions(d,t)=1;
            end
        end
    end
    
    
    if ~mod(t,round(npts/100)), fprintf('  %i%%\n',round(t/npts*100)); end % Print progress
end
toc
fprintf('Number of collisions: %i\n',ncollisions) % Output number of collisions

%% Data analysis and plotting

if doplot==1
    figure('position',[200 200 1600 800])
    radius=tracklength/(2*pi); % Set radius of track
    paint=vectocolor(1:ndrivers)/2; % Choose colors for vehicles
    
    
    %% Plot other quantities

    cksubplot([2 3],[2 1],[70 70],[0 0],[3 -5])
    plot(time,positions'), ylabel('Position (m)'), xlabel('Time (s)')
    cksubplot([2 3],[2 2],[70 70],[0 0],[3 -5])
    plot(time,velocities'), ylabel('Velocity (m/s)'), xlabel('Time (s)')
    cksubplot([2 3],[2 3],[70 70],[0 0],[3 -5])
    plot(time,percheads'), ylabel('Headway (m)'), xlabel('Time (s)')
    
    %% Plot track
    
    h=cksubplot([2 1],[1 1],[80 70],[0 0],[5 -10]); hold on
    grass=[0.2 0.9 0.4]; % Color of grass
    bitumen=0.7*[1 1 1]; % Color of track
    trackwidth=2.5; % Width of track
    axislims=1.2*radius*[-1 1];
    fill(axislims(1)*[-1 1 1 -1],axislims(1)*[-1 -1 1 1],grass)
    phi = (0:0.01:1)'*2*pi;
    xouter = (radius+trackwidth)*sin(phi);
    youter = (radius+trackwidth)*cos(phi);
    fill(xouter,youter,bitumen)
    xinner = (radius-trackwidth)*sin(phi);
    yinner = (radius-trackwidth)*cos(phi);
    fill(xinner,yinner,grass)
    xmid = (radius)*sin(phi);
    ymid = (radius)*cos(phi);
    plot(xmid,ymid,'w--')
    
    %% Plot cars
    
    for t=1:round(0.1/dt):npts
        set(gcf,'currentaxes',h)
        current=positions(:,t);
        
        x=radius*cos(2*pi*current/tracklength);
        y=radius*sin(2*pi*current/tracklength);
        alertness=logical(alert(:,t)); % Show alertness
        w=scatter(x(alertness),y(alertness),100,paint(alertness,:),'filled','y');
        q=scatter(x,y,25,paint,'filled');
        boom=logical(collisions(:,t)); % Show collisions
        z=scatter(x(boom),y(boom),100,'r','filled');
        axis equal; box on
        xlabel('Position (m)')
        ylabel('Position (m)')
        title(sprintf('Current time: %0.1f s',t*dt));
        xlim(axislims)
        ylim(axislims)
        drawnow
        if t<npts, delete(q), delete(w), delete(z), end % Remove points to update
    end
end


disp('Done.')