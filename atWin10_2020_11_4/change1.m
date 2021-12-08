%% Bellhop

filename = 'Sandbottom_2021_6';

% Basic parameters
%c0 = 1500;
%Freq = 1800;

% Environmental Parameters

%water_depth = 1200;
%zmax = 1200;

%rho1 = 1;
%rho2 = 2.1;   %修改参数1 范围先设1.5-2.5 变化间隔0.01  海底密度

%attn1 = 0;
%attn2 = 1; % 修改参数2 范围0.1 - 1 衰减系数

%cb = 1830;  %修改参数3 范围1600 - 2000 间隔10  海底声速


%attn2 = 10;
%theta_0 = 0;
%rho2 = 1.78;
%cb = 1680;
for  rho2 = 1.5 :0.01 : 2.5
    for cb = 1600 : 10 : 2000
        for attn2 = 0 : 0.01 : 1
            for theta_0 = -30 : 1 : 30
                c0 = 1500;
                Freq = 1800;
                water_depth = 1200;
                zmax = 1200;
                rho1 = 1;
                attn1 = 0;
                load('ssp_xctd.mat');
                ssps.svpz = 0 : 5 : water_depth;
                ssps.svp = interp1(ssp_xctd(:,1), ssp_xctd(:,5), ssps.svpz, 'linear', 'extrap');
                globalopt = 'SVW';  % SSP option  N代表SSP插值方式 V代表顶部空间形式 W代表衰减系数单位 W：db/波长
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
                %theta_0 = 0;  % 修改参数4 变化范围 -30 - 30  声源入射角
                theta_range = 60; 
                alpha = [theta_0-theta_range/2,theta_0+theta_range/2];  % 声源出射角
                write_bellhop_env([filename, '.env'], Freq, globalopt, ...
                    media, ...
                    bottom, ...
                    rmax/1e3,nsd,src_depth,nrd,rcv_depth,Pos,Option,Nbeams,alpha);
                if length(Option)>2
                    theta = -90:0.5:90;
                    d = delta_z;
                    f = Freq;
                    c = ssps.svp(1);
                    lambda = c/f;
                    Bp = abs((sin(N*pi*d.*(sind(theta)-sind(theta_0))/lambda))./(sin(pi*d.*(sind(theta)-sind(theta_0))/lambda)));
                    Bpdb= 20*log10(Bp/N);
                    Bpdb(isnan(Bpdb))=0;
                    sbpfil= [ filename '.sbp' ];
                    fid = fopen( sbpfil, 'wt' );
                    fprintf( fid, '%d \n',length(Bpdb) );
                    for ii = 1 : length(Bpdb)
                        fprintf( fid, '%6.2f %6.2f  \t  \n', ...
                            [ theta( ii ) ...
                            Bpdb( ii )  ] );
                    end
                end
                fclose(fid);
                load('scs_bathy_2021.mat'); %海域地形图
                load('GPS.mat')
                [param.src_utmx, param.src_utmy] = deg2utm(Gps.src_lat, Gps.src_lon);
                [param.rcv_utmx, param.rcv_utmy] = deg2utm(Gps.rcv_lat, Gps.rcv_lon);
                src_op = 9190; % 对于实测数据 9181-9271.mat
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
                wd = interp2(xgrid,ygrid,-zz,xv,yv);%wd 正负10% 步长1%
                bathy  = [rv(:)/1e3,wd(:)];
                write_bathy(filename,'C',bathy) % C为地形插值方式
                tic
                bellhop(filename);
                [p_s,rt,zt] = xb_plotshd([filename,'.shd']);%caxis([-100,-50])
                % hold on
                rcv_depth = 51;
                zave = [rcv_depth-5,rcv_depth+5];
                zmin = min(zave(1),wd-1);
                zmax = min(zave(2),wd-1);
                fish_mask = bsxfun(@ge,zt,zmin)&bsxfun(@le,zt,zmax);
                intMid = sum(abs(p_s).*fish_mask,1)./sum(fish_mask,1); % 20*log10(intMid)随距离r变化曲线就是模型输出S
                load('data1_9181_9271.mat')
                tl_interp = interp1(rv,20*log10(intMid),range);
                cp = mse(tl_interp-(p_mf-200));
               % fileID1 = fopen('msewd.txt','a+');
            %    fileID2 = fopen('cpwd1.txt','a+');
                %fprintf(fileID1,'rho2=%.2f,cb=%d,attn2=%.2f,theta_0=%d,cp=%f\n',rho2,cb,attn2,theta_0,cp);
               % fprintf(fileID1,'i=%.2f,cp=%f\n',i,cp);
              %  fprintf(fileID2,'%f\n',cp);
            %    fclose(fileID1);
             %   fclose(fileID2);
            end
        end
    end
end

