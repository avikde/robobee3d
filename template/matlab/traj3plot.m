function traj3plot(ax, ax2, ax3, tt, q, ss, tmax, col)

xyz = q(:,1:3);
ii = find(tt < tmax & tt > 0.01);
plot3(ax, xyz(ii,1), xyz(ii,2), xyz(ii,3), 'color', col)
hold(ax, 'on')
% arrows
jj = ii(1:500:end);

scale = 0; % disable auto scaling
ss = 0.01 * ss; % manual scaling
quiver3(ax, xyz(jj,1), xyz(jj,2), xyz(jj,3), ss(jj,1), ss(jj,2), ss(jj,3), scale, 'color', col)


% traj
plot(ax2, tt(ii), q(ii,4))
hold(ax2, 'on')
plot(ax3, tt(ii), q(ii,5))
hold(ax3, 'on')

end
