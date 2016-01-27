%% Perform computation
r = (0:1:70) *pi/180;
e = (-30:1:30) *pi/180; % epsilon

Naz = (0:10:180)*pi/180;
Nel = (repmat(0,1,numel(Naz)))*pi/180;
[x,z,y] = sph2cart(Naz, Nel, 1);
N = [x; y; z];

theta = 0;
N(3,:) = N(3,:) - 1;
N = [1, 0, 0; 0, cosd(theta), -sind(theta); 0 sind(theta), cosd(theta)]*N;
N(3,:) = N(3,:) + 1;

angles = zeros(numel(r), numel(e), size(N,2));
rcs = zeros(size(angles));

for i = 1:size(angles, 1)
    for j = 1:size(angles, 2)
        for k = 1:size(angles, 3)
            [angle, rc] = q3sim(r(i), N(:,k), e(j));
            assert(isreal(angle))
            angles(i,j,k) = angle;
            rcs(i,j,k) = rc;
        end
    end
end

%% Plot
angles2 = angles;
%angles2(angles2>40) = NaN;
for n = 1:numel(Naz);
    [X, Y] = meshgrid(e, r);
    colormap(flipud(jet(256)));
    h = surf(X*180/pi, Y*180/pi, -angles2(:,:,n), rcs(:,:,n), ...
            'LineStyle', 'none', ...
            'EdgeAlpha', 0.3);
    caxis([0, 1]);
    zlim([-40, 0]);
    
    numberOfXTicks = 5;
    xData = get(h,'XData');
    set(gca,'Xtick',linspace(xData(1),xData(end),numberOfXTicks))
    yData = get(h,'YData');
    set(gca,'Ytick',linspace(yData(1),yData(end),numberOfXTicks))
    set(gca,'fontsize', 16);
    
    xlabel('\epsilon'); ylabel('r'); zlabel('Angular error');
    c = colorbar();
    ylabel(c, 'Reciprocal Condition Number');
    title(['\theta_{azimuth}' sprintf(' = %d', round(Naz(n)*180/pi)) '\circ']);
    set(gcf, 'Position', [676   504   635   434]);
    
    view(-33, 50);
    set(gca, 'ZTickLabel', num2str(-str2double(get(gca,'ZTickLabel'))));
    export_fig(sprintf('q3_erac_Naz%03d_1.pdf', round(Naz(n)*180/pi)), '-transparent');
    view(-160, 50);
    export_fig(sprintf('q3_erac_Naz%03d_2.pdf', round(Naz(n)*180/pi)), '-transparent');
    %break
end
%%
angles2 = angles;
%angles2(angles2>40) = NaN;
for rv = 1:5:numel(r)
    [X, Y] = meshgrid(Naz, e);
    colormap(flipud(jet(256)));
    h = surf(X*180/pi, Y*180/pi, -squeeze(angles2(rv,:,:)), squeeze(rcs(rv,:,:)), ...
             'LineStyle', '-', ...
             'EdgeAlpha', 0.3);
    caxis([0, 1]);
    zlim([-40, 0]);
    
    numberOfXTicks = 5;
    xData = get(h,'XData');
    set(gca,'Xtick',linspace(xData(1),xData(end),numberOfXTicks))
    yData = get(h,'YData');
    set(gca,'Ytick',linspace(yData(1),yData(end),numberOfXTicks))
    set(gca,'fontsize', 16);
    
    ylabel('\epsilon'); xlabel('\theta_{azimuth}'); zlabel('Angular error');
    c = colorbar();
    ylabel(c, 'Reciprocal Condition Number');
    title([sprintf('r = %d', round(r(rv)*180/pi)) '\circ']);
    set(gcf, 'Position', [676   504   635   434]);
    view(-55, 60);
    set(gca, 'ZTickLabel', num2str(-str2double(get(gca,'ZTickLabel'))));
    export_fig(sprintf('q3_neac_r%03d_1.pdf', round(r(rv)*180/pi)), '-transparent');
    view(-160, 60);
    export_fig(sprintf('q3_neac_r%03d_2.pdf', round(r(rv)*180/pi)), '-transparent');
end