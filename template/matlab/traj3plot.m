function traj3plot(tt, xyz, tmax)

ii = find(tt < tmax & tt > 0.01);
plot3(xyz(ii,1), xyz(ii,2), xyz(ii,3))

end
