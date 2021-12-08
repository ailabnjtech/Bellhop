function [p_all,rt,zt] = xb_plotshd(filename)
% plot pressure field in dB 

itheta = 1;
[ PlotTitle, ~, freqVec, ~, ~, Pos, pressure ] = read_shd( filename );
freq = freqVec( 1 );
zt       = reshape(Pos.r.z,[],1);
rt       = reshape(Pos.r.r,1,[]);
% Initialize
p_all = 0;

% 
for isz = 1:size(pressure,2)
    p_one = squeeze( pressure( itheta, isz, :, : ) );
    p_all = p_all+p_one;
end

p_all = p_all/size(pressure,2);
tlt = double(abs(p_all));
tlt( isnan( tlt ) ) = 1e-10;   % remove NaNs
tlt( isinf( tlt ) ) = 1e-10;   % remove infinities

% icount = find( tlt > 1e-37 );        % for stats, only these values count
tlt( tlt < 1e-37 ) = 1e-37;          % remove zeros
tlt_db = 20.0 * log10( tlt );          % so there's no error when we take the log

% tlmed = median( tlt_db( icount ) );    % median value
% tlstd = std( tlt_db( icount ) );       % standard deviation
% tlmax = tlmed + 0.75 * tlstd;       % max for colorbar
% tlmax = 10 * round( tlmax / 10 );   % make sure the limits are round numbers
% tlmin = tlmax - 50;                 % min for colorbar
%% Plot

%rmax = max(Pos.r.r);
%zmax = max(Pos.r.z);
%imagesc(rt,zt,tlt_db)
%colormap( 'jet' )
%caxis( [ -80,-30 ] )
%set( gca, 'xtick',0:5e3:rmax,'xticklabel',0:5:rmax/1e3,'YDir', 'Reverse','Fontsize',16);
%xlim([0,rmax]);ylim([0,zmax])
%xlabel( 'Range(m)' )
%ylabel( 'Depth (m)' );
%set(gcf,'color','w');
%colorbar;

end