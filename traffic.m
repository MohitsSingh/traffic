%% TRAFFIC
% Our totally awesome traffic model.
% Version: 2012nov03

%% Set parameters
tic
rng(239487)
disp('Setting parameters...')
ndrivers=20; % Number of drivers
tracklength=100; % Track length in meters
maxtime=200; % Number of seconds to simulate
dt=0.05; % Timestep size in seconds
npts=maxtime/dt;
time=dt:dt:maxtime;
defaultv=30; % Default velocity in m/s
tau=round(0.5/dt); % Reaction time
C(1)=16.8; % 
C(2)=0.086;
C(3)=C(2)*-25;
C(4)=C(1)*0.913;
alpha=0.5;
randvel=1;
randpos=1;
doplot=0;

vtilda=@(headway) max(0,C(1)*tanh(C(2)*headway+C(3))+C(4));

dx=0.1;
vtildatable=alpha*dt*vtilda(0:dx:200);


%% Set up initial conditions
drivers=struct;
for d=1:ndrivers
    drivers(d).x=(d/ndrivers)*tracklength+randpos*randn; % Initial position
    drivers(d).v=defaultv+randvel*randn; % Initial velocity
    drivers(d).a=0; % Initial acceleration
    drivers(d).state=0; % Alert vs. fucked
    drivers(d).headway=Inf; % Amount of distance to car in front
    drivers(d).infront=mod(d,ndrivers)+1; % Index of driver in front
end


%% Run simulation
disp('Running simulation...')
matrix=zeros(ndrivers,npts);
positions=matrix;
velocities=matrix;
headways=matrix;
states=matrix;
for t=1:npts
    for d=1:ndrivers
        positions(d,t)=drivers(d).x; % Save current position
        velocities(d,t)=drivers(d).v; % Save current position
        headways(d,t)=drivers(d).headway; % Save current headway
        states(d,t)=drivers(d).state; % Save alertness
        
        drivers(d).state=0; % Update state
        if t>tau % Update acceleration
            drivers(d).headway=mod(diff(positions([d drivers(d).infront],t-tau)),tracklength);
            drivers(d).v=(1-alpha*dt)*drivers(d).v+vtildatable(round(drivers(d).headway/dx)+1);
        end
        
        drivers(d).v=drivers(d).v+drivers(d).a; % Update velocity
        drivers(d).x=mod(drivers(d).x+drivers(d).v*dt,tracklength); % Update position
    end
    if ~mod(t,round(npts/100)), fprintf('  %i%%\n',round(t/npts*100)); end
end


if doplot==1
    figure('position',[200 200 1600 800])
    radius=tracklength/(2*pi);
    toc
    
    %% Plot other quantities
    cksubplot([2 3],[2 1],[70 70],[0 0],[3 -5])
    plot(time,positions'), ylabel('Position (m)'), xlabel('Time (s)')
    cksubplot([2 3],[2 2],[70 70],[0 0],[3 -5])
    plot(time,velocities'), ylabel('Velocity (m/s)'), xlabel('Time (s)')
    cksubplot([2 3],[2 3],[70 70],[0 0],[3 -5])
    plot(time,headways'), ylabel('Headway (m)'), xlabel('Time (s)')
    
    %% Plot track
    cksubplot([2 1],[1 1],[80 70],[0 0],[5 -10]) ; hold on
    grass=[0.2 0.9 0.4];
    bitumen=0.7*[1 1 1];
    trackwidth=2.5; % Width of track
    axislims=1.1*radius*[-1 1];
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
        current=positions(:,t);
        x=radius*cos(2*pi*current/tracklength);
        y=radius*sin(2*pi*current/tracklength);
        q=scatter(x,y,'filled','k');
        axis equal; box on
        xlabel('Position (m)')
        ylabel('Position (m)')
        title(sprintf('Current time: %0.1f s',t*dt));
        xlim(axislims)
        ylim(axislims)
        drawnow
        pause(dt/2)
        if t<npts, delete(q), end
    end
end




disp('Done.')