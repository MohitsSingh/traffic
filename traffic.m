%% TRAFFIC
% Our totally awesome traffic model.
% Version: 2012nov03

%% Set parameters
tic
disp('Setting parameters...')
ndrivers=100; % Number of drivers
tracklength=1000; % Track length in meters
maxtime=20; % Number of seconds to simulate
dt=0.02; % Timestep size in seconds
npts=maxtime/dt;
defaultv=30; % Default velocity in m/s

drivers=struct;
for d=1:ndrivers
    drivers(d).x=(d/ndrivers)*tracklength;
    drivers(d).v=defaultv;
end



%% Run simulation
positions=zeros(ndrivers,npts);
for t=1:npts
    for d=1:ndrivers
        positions(d,t)=drivers(d).x;
        drivers(d).x=mod(drivers(d).x+drivers(d).v*dt,tracklength);
    end
end

figure('position',[200 200 800 800]); hold on
radius=tracklength/(2*pi);
toc

%% Plot track
trackwidth=10; % Width of track
phi = (0:0.01:1)'*2*pi;
xouter = (radius+trackwidth)*sin(phi);
youter = (radius+trackwidth)*cos(phi);
fill(xouter,youter,0.7*[1 1 1])
xinner = (radius-trackwidth)*sin(phi);
yinner = (radius-trackwidth)*cos(phi);
fill(xinner,yinner,'w')



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
    xlim(1.1*radius*[-1 1])
    ylim(1.1*radius*[-1 1])
    drawnow
    pause(dt)
    delete(q)
end


disp('Done.')