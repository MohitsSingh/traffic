%% TRAFFIC
% Our totally awesome traffic model.
% Version: 2012nov04

%% Parameters and settings

tic
rng(239487) % Random number seed
disp('Setting parameters...')
ndrivers=6; % Number of drivers
tracklength=300; % Track length in meters
maxtime=40; % Number of seconds to simulate
dt=0.1; % Timestep size in seconds
npts=maxtime/dt;
time=dt:dt:maxtime;
defaultv=20; % Default velocity in m/s
tau=round(0.5/dt); % Reaction time (in number of timesteps)
C(1)=16.8;
C(2)=0.086;
C(3)=C(2)*-25;
C(4)=C(1)*0.913;
alpha=0.5; % Responsivity of drivers
randvel=10; % Randomness in initial velocity in m/s
randpos=20; % Randomness in initial positions in m

psleep = 0.1; % Probability of losing alertness
pwake =  0.1; % Probability of regaining alertness

doplot=1; % Toggle whether to plot
ncollisions=0; % Initialize number of collisions

%% Define Optimal Velocity function, vtilda
vtilda=@(perceivedhead) max(0,C(1)*tanh(C(2)*perceivedhead+C(3))+C(4)); % Taken from paper

dx=0.1; % Spatial step for calculating look-up table
vtildatable=alpha*dt*vtilda(0:dx:tracklength); % Create look-up table for vtilda

%% Set up initial conditions

drivers=struct;
for d=1:ndrivers
    drivers(d).x=(d/ndrivers)*tracklength+randpos*randn; % Initial position
    drivers(d).v=defaultv+randvel*randn; % Initial velocity
    drivers(d).state=0; % Alert vs. fucked
    drivers(d).realhead=tracklength/ndrivers; % Amount of distance to car in front
    drivers(d).perceivedhead=tracklength/ndrivers; % Amount of distance to car in front
    drivers(d).infront=mod(d,ndrivers)+1; % Index of driver in front
    drivers(d).behind=mod(d-2,ndrivers)+1; % Index of driver behind
end

%% Run simulation

% Initialize data matrices
disp('Running simulation...')
matrix=zeros(ndrivers,npts);
positions=matrix;
velocities=matrix;
realheads=matrix;
perceivedheads=matrix;
states=matrix;
collisions=matrix;

% Solve difference equations across all time points
for t=1:npts
    
    
    % Solve motion of each driver at timepoint 'd'
    for d=1:ndrivers
        if t>tau % Update velocity
            sr = rand(1); % Random variable for determing transitions of Markov process
            if drivers(d).state==0
                if sr<psleep,
                    drivers(d).state=1;
                else
                    drivers(d).state=0;
                end
            else
                if sr<pwake,
                    drivers(d).state=0;
                else
                    drivers(d).state=1;
                end
            end
            
            %drivers(d).state=0; % Update state
            if drivers(d).state == 0,
                drivers(d).perceivedhead=mod(diff(positions([d drivers(d).infront],t-tau)),tracklength); % Update perceived headway if driver is alert
            end
            drivers(d).v=(1-alpha*dt)*drivers(d).v+vtildatable(round(drivers(d).perceivedhead/dx)+1); % Update velocity as a function of perceived headway
        end
        drivers(d).x=mod(drivers(d).x+drivers(d).v*dt,tracklength); % Update position
        
        positions(d,t)=drivers(d).x; % Save current position
        velocities(d,t)=drivers(d).v; % Save current position
        realheads(d,t)=drivers(d).realhead; % Save current actual headway
        perceivedheads(d,t)=drivers(d).perceivedhead; % Save current perceive dhead
        states(d,t)=drivers(d).state; % Save state
    end
    
    % Use actual headway to check for drivers overtaking one another; treat these events as collisions
    for d=1:ndrivers
        % Update positions
        if t>1
  
            drivers(d).realhead=mod(diff(positions([d drivers(d).infront],t-1)),tracklength);
            % Check if collision occurred. If it did, reorder drivers and
            % alter velocities of collided vehicles
            count=0;
            while drivers(d).realhead>tracklength-2*dt*(C(1)+C(4))
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
    plot(time,perceivedheads'), ylabel('Headway (m)'), xlabel('Time (s)')
    
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
    
    for t=1:npts
        set(gcf,'currentaxes',h)
        current=positions(:,t);
        
        x=radius*cos(2*pi*current/tracklength);
        y=radius*sin(2*pi*current/tracklength);
        alertness=logical(states(:,t)); % Show alertness
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
%         pause(dt) % Choose rate of animation
        if t<npts, delete(q), delete(w), delete(z), end % Remove points to update
    end
end


disp('Done.')