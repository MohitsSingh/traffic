%% TRAFFIC
% Our totally awesome traffic model.
% Version: 2012nov03

%% Set parameters
tic
disp('Setting parameters...')
ndrivers=20; % Number of drivers
tracklength=1000; % Track length in meters
maxtime=20; % Number of seconds to simulate
dt=0.1; % Timestep size in seconds
npts=maxtime/dt;
defaultv=30; % Default velocity in m/s
tau=round(0.3/dt); % Reaction time
C(1)=16.8; % 
C(2)=0.086;
C(3)=C(2)*-25;
C(4)=0.913;
alpha=0.5;

vtilda=@(headway) max(0,C(1)*tanh(C(2)*headway+C(3))+C(4));


%% Set up initial conditions
drivers=struct;
for d=1:ndrivers
    drivers(d).x=(d/ndrivers)*tracklength; % Initial position
    drivers(d).v=0;% defaultv+5*randn; % Initial velocity
    drivers(d).a=0; % Initial acceleration
    drivers(d).state=0; % Alert vs. fucked
    drivers(d).headway=Inf; % Amount of distance to car in front
    drivers(d).infront=mod(d,ndrivers)+1; % Index of driver in front
end


%% Run simulation
positions=zeros(ndrivers,npts);
velocities=zeros(ndrivers,npts);
for t=1:npts
    for d=1:ndrivers
        positions(d,t)=drivers(d).x; % Save current position
        velocities(d,t)=drivers(d).v; % Save current position
        
        drivers(d).state=0; % Update state
        if t>tau % Update acceleration
            drivers(d).headway=diff(positions([d drivers(d).infront],t-tau));
            drivers(d).v=drivers(d).v+alpha*dt*(-drivers(d).v+vtilda(drivers(d).headway));
        end
        
        drivers(d).v=drivers(d).v+drivers(d).a; % Update velocity
        drivers(d).x=mod(drivers(d).x+drivers(d).v*dt,tracklength); % Update position
        
    end
end

figure('position',[200 200 800 800]); hold on
radius=tracklength/(2*pi);
toc


%% Plot track
grass=[0.2 0.9 0.4];
bitumen=0.8*[1 1 1];
trackwidth=10; % Width of track
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
%     pause(dt)
    if t<npts, delete(q), end
end


disp('Done.')