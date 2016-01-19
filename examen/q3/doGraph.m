%% Perform computation
r = (0:1:80) *pi/180;
e = (0:1:30) *pi/180; % epsilon

Naz = (0:10:180)*pi/180;
Nel = (repmat(0,1,numel(Naz)))*pi/180;
[x,z,y] = sph2cart(Naz, Nel, 1);
N = [x; y; z];

angles = zeros(numel(r), numel(e), size(N,2));
rcs = zeros(size(angles));

for i = 1:size(angles, 1)
    for j = 1:size(angles, 2)
        for k = 1:size(angles, 3)
            [angle, rc] = sim(r(i), N(:,k), e(j));
            angles(i,j,k) = angle;
            rcs(i,j,k) = rc;
        end
    end
end

%% Plot
angles2 = angles;
angles2(angles2>40) = NaN;
for n = 1:numel(Naz);

    [X, Y] = meshgrid(e, r);
    colormap(flipud(jet(256)));
    surfc(X*180/pi, Y*180/pi, angles2(:,:,n), rcs(:,:,n));
    caxis([0, 1]);
    zlim([0, 40]);
    xlabel('epsilon'); ylabel('r'); zlabel('Angular error');
    c = colorbar();
    ylabel(c, 'Reciprocal Condition Number');
    view(-100, 30);

    % TODO: change to pdf!
    export_fig(sprintf('q3_erac_Naz%03d.png', round(Naz(n)*180/pi)), '-transparent', '-m2');
end
%%
angles2 = angles;
angles2(angles2>40) = NaN;
for rv = 1:numel(r)
    [X, Y] = meshgrid(Naz, e);
    colormap(flipud(jet(256)));
    surfc(X*180/pi, Y*180/pi, squeeze(angles2(rv,:,:)), squeeze(rcs(rv,:,:)));
    caxis([0, 1]);
    zlim([0, 40]);
    ylabel('epsilon'); xlabel('theta azimuth'); zlabel('Angular error');
    c = colorbar();
    ylabel(c, 'Reciprocal Condition Number');
    view(-15, 30);
    %break
    
    export_fig(sprintf('q3_neac_r%03d.png', round(r(rv)*180/pi)), '-transparent', '-m2');
    break
end