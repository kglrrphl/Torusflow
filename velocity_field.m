%% velocity field

path = '';
flag1 = 1; % 1 = create new coeffs
flag2 = 0; % 1 = write files
flag3 = 0; % 1 = save picture
calclim = [0, 2*pi]; % calculation range

%% set up coeffs, stream function, velocity field
rho = (sqrt(5)-1)/2;
f = 1.5;
num = 1; % choose from below

if (flag1 == 0)
    load([path,'parameters.mat'])
else  
    switch num
        case 1 %symmetric
            dim = 2;
            k = 0:dim - 1;
            cxcy = [0,1/sqrt(1 + 1); -1/sqrt(rho^2 + 1),0];
            amp = sqrt(rho^2 + 1);
            cxcy = cxcy*amp;
            phix = [pi/2, 0];
            phiy = [pi/2, 0];

        case 2 %asymmetric 
            dim = 2;
            k = 0:dim - 1;
            cxcy = [0,1;-1,-1];
            cxcy = cxcy/1.5843;
            phix = [pi/2, pi/2];
            phiy = [pi/2, -pi/2];
            %phix = [pi/2, 0];
            %phiy = [pi/2, -1];

        case 3 % random
            dim = 15;
            k = 0:dim - 1;
            ks = sqrt(k);
            kk = kron(1./ks', 1./ks); kk(:, 1) = 1./ks'; kk(1, :) = 1./ks; kk(1, 1) = 1; 
            cxcy = reshape(randcoeff(dim^2, 1), dim, dim);
            amp = sqrt(sum(sum(abs(cxcy)*k'))*sum(sum(abs(cxcy')*k')))/sqrt(rho);
            cxcy = cxcy/amp;
            phix = [pi/2, randcoeff(dim - 1, 2*pi)];
            phiy = [pi/2, randcoeff(dim - 1, 2*pi)];
            
        case 4 %perturbation to asymmetric
            dim = 20;
            k = 0:dim - 1;
            cxcy = zeros(dim, dim);
            cxcy(1:2, 1:2) = [0,1;1,1];
            cxcy(end, end) = 0.02;
            cxcy = cxcy/sqrt(2.615);
            phix = [pi/2, zeros(1, dim - 1)];
            phiy = [pi/2, zeros(1, dim - 1)];
            
        case 5 %perturbation to special
            dim = 20;
            k = 0:dim - 1;
            cxcy = zeros(dim, dim);
            cxcy(1:2, 1:2) = [0,1/sqrt(2); -1/sqrt(rho^2+1),0];
            amp = sqrt(rho^2+1);
            cxcy = cxcy*amp;
            cxcy(end, end) = 0.003;
            phix = [pi/2, zeros(1, dim - 1)];
            phiy = [pi/2, zeros(1, dim - 1)];
            
        case 6 %symmetric extension 1
            dim = 3;
            k = 0:dim - 1;
            cxcy = [0, 0.8, 0.4; -1, 0, 0; 0, 0, 0];
            phix = [pi/2, 0, 0];
            phiy = [pi/2, 0, 0];
            
        case 7 %symmetric extension 2
            dim = 3;
            k = 0:dim - 1;
            cxcy = [0, 0.8, 0.4; -1, 0, 0; -0.5, 0, 0];
            cxcy = cxcy*8/9;
            phix = [pi/2, 0, 0];
            phiy = [pi/2, 0, 0];
            
        case 8 %simpler asymmetric
            dim = 2;
            k = 0:dim - 1;
            cxcy = [0,0;0,1/rho];
            phix = [pi/2, 0];
            phiy = [pi/2, 0];
            
        case 9 %asymmetric extension
            dim = 3;
            k = 0:dim - 1;
            cxcy = [0, 0, 0; 0, 1, 0; 0, 0, 0.5];
            cxcy = cxcy*(1+rho)/2;
            phix = [pi/2, pi, 0];
            phiy = [pi/2, pi/2, pi/2];
            
        case 10 %mixed
            dim = 3;
            k = 0:dim - 1;
            cxcy = [0, 0.8, 0.4; -1, 1, 0; -0.5, 0, 0];
            cxcy = cxcy*0.667;
            phix = [pi/2, 0.6, 1.2];
            phiy = [pi/2, 1, 2];    
            
        case 11 %sawtooth
            dim = 13;
            k = 0:dim - 1;
            saw = @(n) 2/pi*(-1).^(n+1)./n;
            cxcy = zeros(dim,dim); cxcy(1,2:end) = saw(1:dim-1); cxcy(2:end,1) = saw(1:dim-1);
            phix = [pi/2, zeros(1,dim-1)];
            phiy = [pi/2, zeros(1,dim-1)];
            
        case 12 %sawtooth
            dim = 6;
            k = 0:dim - 1;
            saw = @(n) 2/pi*(-1).^(n)./n.^2;
            cxcy = zeros(dim,dim); cxcy(1,2:end) = saw(1:dim-1); cxcy(2:end,1) = saw(1:dim-1);
            phix = [pi/2, pi/2*ones(1,dim-1)];
            phiy = [pi/2, pi/2*ones(1,dim-1)];
            
        case 13 %sawtooth
            dim = 5;
            k = 0:dim - 1;
            saw = @(n,c) 2/(rho*pi)./(n.^3.*sqrt(c^2+n.^2)); phi = @(n,c) atan(c./n);
            cxcy = zeros(dim,dim); cxcy(1,2:end) = saw(1:dim-1,1); cxcy(2:end,1) = -saw(1:dim-1,rho);
            phix = [pi/2, phi(1:dim-1,rho)];
            phiy = [pi/2, phi(1:dim-1,1)];
            
                        
        case 14 %square
            dim = 8;
            k = 0:dim - 1;
            saw = @(n,c) 2/(rho*pi)./(n.^3.*sqrt(c^2+n.^2)); phi = @(n,c) atan(c./n);
            cxcy = zeros(dim,dim); cxcy(1,2:2:end) = saw(1:2:dim-1,1); cxcy(2:2:end,1) = -saw(1:2:dim-1,rho);
            phix = [pi/2, phi(1:dim-1,rho)];
            phiy = [pi/2, phi(1:dim-1,1)];
    end
end

psi = @(x,y) (-x + rho*y) + f*sin(k*x + phix)*cxcy*sin(k*y + phiy)';
v = @(r) [rho + f*(sin(k*r(1) + phix)*cxcy*(cos(k*r(2) + phiy).*k)');
         1 - f*(sin(k*r(2) + phiy)*cxcy'*(cos(k*r(1) + phix).*k)')];

     
%% save coeffs
if (flag2 == 1)
    func = @(x) num2str(x,'%7f');
    txt = ['double cxcy[',num2str(dim),'][',num2str(dim),'] = {',...
           strjoin(arrayfun(func,cxcy','UniformOutput',false),','),'};\n',...
           'vector<double> phix = {',strjoin(arrayfun(func,phix,'UniformOutput',false),','),'};\n',...
           'vector<double> phiy = {',strjoin(arrayfun(func,phiy,'UniformOutput',false),','),'};'];

    fid = fopen([path,'coefficients.txt'],'wt');
    fprintf(fid, txt);
    fclose(fid);

    save([path,'parameters.mat'],'cxcy','phix','phiy','k','dim')
    disp('coeffs saved')
end


%% calc stagnation points

pr = calclim(2) - calclim(1);
np = numel(k)*10;
inc = linspace(calclim(1)-1, calclim(2)+1, np);
sol = zeros(0,2);
rot = [];
val = [];
options = optimset('Display','off','TolFun',1e-12,'FinDiffType','central','MaxIter',1e2,'Tolx',1e-8);

for i = 1:np
    for j = 1:np
        
        [root,fval,ef,~,J] = fsolve(v, [inc(i), inc(j)], options);
        
        if (ef >0)
            
            root = map(root, calclim(1), calclim(2));
            root = round(root,3);
                        
            if (sum(all(ismember(root,sol,'rows'), 2)) == 0)
                
                sol(end+1,:) = root;
                val(end+1,:) = fval;
                
                if (abs(imag(eig(J))) < 1e-2)
                    rot(end+1) = 0;
                else
                    aux = v([root(1)-0.05,root(2)]);
                    rot(end+1) = sign(aux(2));
                end    
            end 
        end
    end
end


%% calc separatrices

lvls = [];
stags = sol(rot==0,:);
for i = 1:size(stags,1)
    lvls(end+1) = psi(stags(i,1), stags(i,2));
    if sum(lvls(end) == lvls(1:end-1)) > 0
        lvls(end) = [];
    end
end

% use this if separatrices form neighbouring cells are needed
% aux1 = [0,2*pi,-2*pi,0,   0,    2*pi,-2*pi,2*pi, -2*pi];
% aux2 = [0,0,   0,    2*pi,-2*pi,2*pi,-2*pi,-2*pi,2*pi];
% for i = 1:size(stags,1)
%     for j = 1:9
%         lvls(end+1) = psi(stags(i,1)+aux1(j), stags(i,2)+aux2(j));
%         if sum(lvls(end) == lvls(1:end-1)) > 0
%             lvls(end) = [];
%         end
%     end
% end

lvls = sort(lvls);
if numel(lvls) == 1
    lvls(2) = lvls(1);
end   


%% isolines of psi

vec = linspace(calclim(1), calclim(2), 2000);
[xxx, yyy] = meshgrid(vec);
Z = zeros(size(xxx));

for i = 1:numel(vec)
    for j = 1:numel(vec)
        Z(i, j) = psi(vec(j), vec(i));
    end
end


%% plot streamlines and separatrices

h = figure;
hold on
contour(xxx, yyy, Z, 10, 'linecolor', [1,0,1]*0.5, 'linewidth', 0.5)
contour(xxx, yyy, Z, lvls, 'linecolor', 'k', 'linewidth', 0.5)

for i = 1:numel(rot)
    if (rot(i) == 1)
        mark = 'r>';
    elseif (rot(i) == -1)
        mark = 'r<';
    else
        mark = 'ro';
    end
plot(sol(i,1), sol(i,2), mark, 'markersize',6)
end
hold off

%xlim([-pi, pi]); ylim([-pi,pi])

% title({['Streamlines with f = ', num2str(f)];
%     [num2str(sum(rot>0)), ' clockwise and ', num2str(sum(rot<0)), ' counter-clockwise eddies'];
%     [num2str(sum(rot==0)), ' saddle points']})

box on
%xlabel('x'); ylabel('y')
set(gca,'xtick',[],'ytick',[])
%xticklabels({'0', '2\pi'}); yticklabels({'2\pi'})
pbaspect([1 1 1])
if (flag3 == 1)
    set(h, 'Color', 'w');
    export_fig([path,'field.pdf'], h)
end