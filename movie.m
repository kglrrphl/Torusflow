path = '';
rng = 6*pi; %window size (1 cell = 2*pi)


%% load values
x = load([path, 'vid_x.txt']);
y = load([path, 'vid_y.txt']);
t = load([path, 'vid_t.txt']);

x = x + 2*pi*load([path, 'vid_ox.txt']);
y = y + 2*pi*load([path, 'vid_oy.txt']);
xmean = mean(x,1); ymean = mean(y,1);

x = map(x, 0, rng); y = map(y, 0, rng);
xmean = map(xmean, 0, rng); ymean = map(ymean, 0, rng);

nfr = size(x,2);
im = cell(nfr);


%% calc isolines
vec = linspace(0, rng, 300);
[xxx, yyy] = meshgrid(vec);
Z = zeros(size(xxx));
for i = 1:numel(vec)
    for j = 1:numel(vec)
        Z(i, j) = psi(vec(j), vec(i));
    end
end


%% gather frames
co = [0,0,1]*0.8;
fig = figure;
set(fig, 'Color', 'w')
set(fig, 'Position', get(0,'Screensize'));
for i = 1:nfr
    contour(xxx, yyy, Z, lvls, 'linecolor', 'k')
    hold on
    plot(x(:, i), y(:, i), '.', 'MarkerEdgeColor', co, 'MarkerFaceColor', co)
    plot(xmean(i), ymean(i), 'o', 'MarkerFaceColor', 'r')
    hold off

    title(['t = ', num2str(t(i), '%.1f'), ' s'])
    xlim([0,rng]); ylim([0,rng])
    set(gca,'xtick',[],'ytick',[])
    pbaspect([1 1 1])
    drawnow
    frame = getframe(fig);
    im{i} = frame2im(frame); 
end
close(fig);


%% create gif
% for i = 1:nfr
%     [A, map] = rgb2ind(im{i},256);
%     if i == 1
%         imwrite(A, map, filename, 'gif', 'LoopCount', Inf, 'DelayTime',3);
%     else
%         imwrite(A, map, filename, 'gif', 'WriteMode', 'append', 'DelayTime', 0.1);
%     end
% end


%% create video
vid = VideoWriter([path, 'movie.avi']);
vid.FrameRate = 30;
vid.Quality = 40;
open(vid);
for i = 1:vid.FrameRate
writeVideo(vid, im{1});
end
for i = 1:nfr
    writeVideo(vid, im{i});
end
close(vid);

