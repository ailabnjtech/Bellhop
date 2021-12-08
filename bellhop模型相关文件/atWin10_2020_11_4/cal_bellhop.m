%% Bellhop

filename = 'Sandbottom_2021_6';

% Basic parameters
c0 = 1500;
Freq = 1800;

% Environmental Parameters

water_depth = 1200;
zmax = 1200;

rho1 = 1;
rho2 = 2;   %参数1

attn1 = 0;
attn2 = 10; % 参数2

cb = 2000;  %参数3

%% SSP

fprintf('>>> Load sound speed .......                  ')

load('I:\SCSEx21\ctd&xctd\xctd\ssp_xctd.mat');
ssps.svpz = 0 : 5 : water_depth;
ssps.svp = interp1(ssp_xctd(:,1), ssp_xctd(:,5), ssps.svpz, 'linear', 'extrap');
fprintf('Done !\n')
figure(3); clf(3);
plot(ssps.svp, ssps.svpz,'linewidth',1.5);
set(gca, 'YDir', 'reverse');
set(gcf,'color','w')
xlabel('Sound Speed (m/s)');
ylabel('Water Depth (m)');

%% write Environmental file

fprintf('>>> Writing Environmental file.......         ')
globalopt = 'SVW';  % SSP option  N代表SSP插值方式 V代表顶部空间形式，V为真空，W代表衰减系数单位 W：db/波长

media{1} = struct('sigma', 0, 'zmax', water_depth, 'svpz', ssps.svpz, ...
    'svp', ssps.svp, 'svs', zeros(size(ssps.svp)), 'rho', rho1 * ones(size(ssps.svp)), ...
    'attnp', attn1 * ones(size(ssps.svp)), 'attns', zeros(size(ssps.svp)));

bottom = struct('opt', 'A*', 'sigma', 0, 'svpz', water_depth, 'svp', cb, ...
    'svs', 0, 'rho', rho2, 'attnp', attn2, 'attns', 0);

rmax = 20e3;  % Rmax 计算的最远距离 单位m 写入环境文件时要取km

N = 10;
zs_min = 48.2;
zs_max = 51.8;
delta_z = (zs_max-zs_min)/(N-1);
zs = zs_min + delta_z*(N-1)/2;
nsd = 1;
src_depth = zs;

dz = c0/Freq/8;
nrd = round(water_depth/dz)+1; % Number of revceiver  与dz对应起来
rcv_depth = [0, water_depth]; % rcv depth 接收深度

dr = 50;
Pos.r.range = 0:dr/1e3:rmax/1e3; % 距离范围 dr = 15m
Pos.r.depth = 0:dz:zmax; % 深度范围
Pos.s.depth = zs;    % 声源深度

Option = 'CG*';

Nbeams = 500;
alpha = [-25,-25];  % 声源出射角 需要修改

write_bellhop_env([filename, '.env'], Freq, globalopt, ...
    media, ...
    bottom, ...
    rmax/1e3,nsd,src_depth,nrd,rcv_depth,Pos,Option,Nbeams,alpha);

fprintf('Done !\n')

%% write source beam pattern file(.sbp)
if length(Option)>2
    theta_0 = (alpha(1)+alpha(2))/2; % theta_0 就是声源出射角度
    theta = -90:0.5:90;
    d = delta_z;
    f = Freq;
    c = ssps.svp(1);
    lambda = c/f;
    Bp = abs((sin(N*pi*d.*(sind(theta)-sind(theta_0))/lambda))./(sin(pi*d.*(sind(theta)-sind(theta_0))/lambda)));
    Bpdb= 20*log10(Bp/N);
    Bpdb((length(Bpdb)+1)/2)=0;
    figure(11);clf(11)
    plot(theta,Bpdb);
    xlabel('\theta');ylabel('source beam pattern')
    set(gcf,'color','w');set(gca,'fontsize',12)
    xlim([-90,90]);
    %     saveas(figure(11),'src_beam_pattern.png')
    
    sbpfil= [ filename '.sbp' ];
    fid = fopen( sbpfil, 'wt' );
    fprintf( fid, '%d \n',length(Bpdb) );
    
    % x = mean(Bpdb(1:181));
    for ii = 1 : length(Bpdb)
        fprintf( fid, '%6.2f %6.2f  \t  \n', ...
            [ theta( ii ) ...
            Bpdb( ii )  ] );
    end
end

%% write bathy file
fprintf('>>> Writing bathy file.......                 ')

load('J:\SCSEx2021\scs_bathy_2021.mat'); %海域地形图
load('J:\SCSEx2021\GPS.mat')
[param.src_utmx, param.src_utmy] = deg2utm(Gps.src_lat, Gps.src_lon);
[param.rcv_utmx, param.rcv_utmy] = deg2utm(Gps.rcv_lat, Gps.rcv_lon);
%
% xmin = min(min(param.src_utmx,param.rcv_utmx))-5e3;
% xmax = max(max(param.src_utmx,param.rcv_utmx))+5e3;
% ymin = min(min(param.src_utmy,param.rcv_utmy))-5e3;
% ymax = max(max(param.src_utmy,param.rcv_utmy))+5e3;
%
% arbox = [xmin,xmax,ymin,ymax];
% inx = xgrid >= arbox(1) & xgrid <= arbox(2);
% iny = ygrid >= arbox(3) & ygrid <= arbox(4);
%
src_op = 9190; % 对于实测数据 9181-9271.mat
% src_op = 9310; % 对于实测数据 9281-9348.mat
% 根据不同的src-rcv坐标 可以得到他们之间的地形
src_x0 = param.src_utmx(src_op);
src_y0 = param.src_utmy(src_op);
rcv_x0 = param.rcv_utmx(src_op);
rcv_y0 = param.rcv_utmy(src_op);

angle = atan2d(rcv_y0-src_y0,rcv_x0-src_x0);

xo = src_x0;
yo = src_y0;
rv = 0 : dr : rmax;
xv = xo + rv * cosd(angle);
yv = yo + rv * sind(angle);
wd = interp2(xgrid,ygrid,-zz,xv,yv);
bathy  = [rv(:)/1e3,wd(:)];

write_bathy(filename,'C',bathy) % C为地形插值方式

fprintf('Done !\n')

%% Calculate G(r|rs, f)

fprintf('>>> Calculate Green function.......           ')

tic
system(sprintf('H:\\atWin10_2020_11_4\\atWin10_2020_11_4\\windows-bin-20201102\\bellhop.exe %s', filename));
fprintf('Done !')
fprintf('      耗时：%.2f min \n',toc/60)

%% Pressure field
fprintf('>>> Plot pressure field.......                ')

figure(5);clf(5);
[p_s,rt,zt] = xb_plotshd([filename,'.shd']);caxis([-100,-50])
hold on
plot(rv,wd,'w-','linewidth',2)
fprintf('Done !\n')

hold on; 
plot(range,p_mf-200,'k*');
