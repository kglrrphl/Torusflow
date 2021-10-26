%% plot from c++ output
flag = 1; % 1 = write files
flagc = 1; % 1 = compute C_int

%% load values
path = '';

zeta = load([path,'variance.txt']);
meanv = load([path,'mean.txt'])/(2*pi);
tau = load([path,'tscale.txt']);
t = load([path,'time.txt']);
p_nc = load([path0,'p_nc.txt']);


%% plot zeta

%fits
ttlin = 1000;
ttlog = 100;

%lin fit
mlin = meanv(t > ttlin & t < 1e6);
zlin = zeta(t > ttlin & t < 1e6);
fitlin = @(a,x) a(2).*x.^a(1);
clin = polyfit(log(mlin),log(zlin),1);
clin(2) = exp(clin(2));

%log fit
mlog = meanv(t > ttlog & t < 1e6);
zlog = zeta(t > ttlog & t < 1e6);
fitlog = @(c,x) c(1).*log(c(2).*x).^2;
clog = real(nlinfit(mlog, zlog, fitlog, [0.5, 1]));


h1 = figure;
hold on
g1 = plot(meanv, zeta);
g2 = plot(mlog, 1.5*fitlog(clog,mlog), '--');
%g3 = plot(mlin, fitlin(clin,mlin), '--');
title(['\zeta for p_{nc} = ', num2str(round(p_nc,3))])
legend(g2,'$\sim \log^2 n$', 'location','southeast', 'interpreter', 'latex')
%legend(g3,'$\sim t^{',num2str(clin(1)),'}$'], 'location','southeast','interpreter','latex')
xlim([meanv(1), meanv(end)]);
ylim([1e-1,1e3]); 
hold off; 
xlabel('mean number of revolutions')
ylabel('\zeta');
set(gca, 'XScale', 'log') 
set(gca, 'YScale', 'log')

grid on

if (flag == 1)
    set(h1, 'Color', 'w');
    export_fig([path,'zeta.pdf'], h1)
end


%% plot correlation function
cc = load([path,'correlation.txt']);
ttt = (1:numel(cc));
ttt = ttt(ttt < t(end)*0.1);
cc = cc(2:numel(ttt)+1);

h2 = figure;
plot(ttt*tau/(p_nc*2*pi), cc);
set(gca, 'XScale', 'log');
hh = refline(0, 0); hh.Color = 'k';
xlim([0.2,ttt(end)*tau/(p_nc*2*pi)]);
xlim([meanv(1), meanv(end)/10])
ylim([-0.5, 0.5])
xlabel('mean number of revolutions'); 
ylabel('C');
pbaspect([3 1 1])

if (flag == 1)
    set(h2, 'Color', 'w');
    export_fig([path,'correlation.pdf'], h2)
end


%% plot C_int
if(flagc ==1)
    tint = unique(floor(logspace(0, log10(ttt(end)), 100)));
    cint = zeros(1,numel(tint));

    for i = 1:numel(cint)
        cint(i) = (0.5+sum(cc(1:tint(i)).^2)-cc(tint(i)).^2/2)/tint(i);
    end

    tfit = tint(tint > 1e2);
    cfit = cint(tint > 1e2);
    d = polyfit(log(tfit), log(cfit), 1);

    h3 = figure;
    plot(tint, cint)
    hold on
    g4 = plot(tfit, 2*exp(d(2))*tfit.^d(1));
    hold off
    legend(g4,['$C_{int} \sim t^{ ', num2str(d(1),2), '}$'], 'location', 'northeast', 'interpreter', 'latex')
    set(gca, 'XScale', 'log', 'YScale', 'log');
    xlabel('t'); ylabel('C_{int}')
end


%% plot power spectrum
ps = load([path,'spectrum.txt']);
fff = pi*(0:numel(ps)-1)./(numel(ps));

h4 = figure;
plot(fff,ps)
xlim([0, 0.3]); %ylim([1e-2, 1e4])
xlabel('\omega'); ylabel('S')
set(gca, 'YScale', 'log');

if (flag == 1)
    set(h4, 'Color', 'w');
    export_fig([path,'spectrum.pdf'], h4)
end
