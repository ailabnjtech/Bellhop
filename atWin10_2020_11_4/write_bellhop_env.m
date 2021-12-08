function write_bellhop_env(FILENAME, FREQ, GLOBALOPT,...
    MEDIA, ...
    ...nmesh1,sigma1,zmax1,svpz1,svp1,svs1,rho1,attnp1,attns1,...
    ...nmesh2,sigma2,zmax2,svpz2,svp2,svs2,rho2,attnp2,attns2,...
    BOTTOM, ...
    ...botopt,botsigma,botsvpz,botsvp,botsvs,botrho,botattnp,botattns,...
    RMAX,nsd,src_depth,nrd,rcv_depth,Pos,Option,Nbeams,alpha)

fid = fopen(FILENAME,'w');

%% Title
context = sprintf('''%s''',FILENAME);
fprintf(fid, '%s\n', mkstr(context, '! TITLE', 57));

%% freq
context = sprintf('%.3f', FREQ);
fprintf(fid, '%s\n', mkstr(context, '! FREQ (Hz)', 57));

%% nmedia
nmedia = length(MEDIA);
context = sprintf('%d', nmedia);
fprintf(fid, '%s\n', mkstr(context, '! NMEDIA', 57));

%% options
context = sprintf('''%s''', GLOBALOPT);
fprintf(fid, '%s\n', mkstr(context, '! OPTIONS', 57));

%% Layers

for ii = 1:length(MEDIA)
    % nmesh sigma zmax
    context = sprintf('0  %.3f  %.3f', MEDIA{ii}.sigma, MEDIA{ii}.zmax);
    fprintf(fid, '%s\n', mkstr(context, '! NMESH  SIGMA  ZMAX', 57));
    
    % z cp cs rho attnp attns
    context = sprintf('%8.3f %8.3f %4.3f %4.3f %4.3f %4.3f', ...
        MEDIA{ii}.svpz(1), MEDIA{ii}.svp(1), MEDIA{ii}.svs(1), MEDIA{ii}.rho(1), MEDIA{ii}.attnp(1), MEDIA{ii}.attns(1));
    fprintf(fid, '%s\n', mkstr(context, '! Z  CP  CS  RHO(g/cm3)  ATTN_P  ATTN_S', 57));
    
    for jj = 2:length(MEDIA{ii}.svpz)
        fprintf(fid, '%8.3f %8.3f /\n', MEDIA{ii}.svpz(jj), MEDIA{ii}.svp(jj));
    end
    
end

% %% Layer 1
% % nmesh, sigma, zmax
% context = sprintf('%.0f  %.1f  %.1f', nmesh1, sigma1, zmax1);
% fprintf(fid, '%s\n', mkstr(context, '! NMESH  SIGMA  ZMAX(NSSP)', 57));
%
% % z cp cs rho ap as
% context = sprintf('%8.1f %8.2f %4.1f %4.1f %4.1f %4.1f', ...
%     svpz1(1), svp1(1), svs1(1), rho1(1), attnp1(1), attns1(1));
% fprintf(fid, '%s\n', mkstr(context, ...
%     '! Z  CP  CS  RHO(g/cm3)  ATTN_P  ATTN_S', 57));
%
% % for ii = 2:length(svpz1)
% %     fprintf(fid,'%8.1f %8.2f %4.1f %4.1f %4.1f %4.1f\n', [svpz1(ii) svp1(ii) svs1(ii) rho1(ii) attnp1(ii) attns1(ii) ]);
% % end
% for ii = 2:length(svpz1)
%     fprintf(fid,'%8.1f %8.2f /\n', [svpz1(ii) svp1(ii)]);
% end
%
% %% Layer 2 ~ N
% % nmesh, sigma, zmax
% context = sprintf('%.0f  %.1f  %.1f', nmesh2, sigma2, zmax2);
% fprintf(fid, '%s\n', mkstr(context, '! NMESH  SIGMA  ZMAX(NSSP)', 57));
%
% % z cp cs rho ap as
% context = sprintf('%8.1f %8.2f %4.1f %4.1f %4.1f %4.1f', ...
%     svpz2(1), svp2(1), svs2(1), rho2(1), attnp2(1), attns2(1));
% fprintf(fid, '%s\n', mkstr(context, ...
%     '! Z  CP  CS  RHO(g/cm3)  ATTN_P  ATTN_S', 57));
%
% % for ii = 2:length(svpz2)
% %     fprintf(fid,'%8.1f %8.2f %4.1f %4.1f %4.1f %4.1f\n', [svpz2(ii) svp2(ii) svs2(ii) rho2(ii) attnp2(ii) attns2(ii) ]);
% % end
% for ii = 2:length(svpz2)
%     fprintf(fid,'%8.1f %8.2f /\n', [svpz2(ii) svp2(ii)]);
% end

%% Bottom
% opt sigma
context = sprintf('''%s''  %.3f', BOTTOM.opt, BOTTOM.sigma);
fprintf(fid, '%s\n', mkstr(context, '! BOTOPT, SIGMA (m)', 57));

% z cp cs rho ap as
context = sprintf('%8.3f %8.3f %4.3f %4.3f %4.3f %4.3f', [BOTTOM.svpz BOTTOM.svp BOTTOM.svs BOTTOM.rho BOTTOM.attnp BOTTOM.attns]);
fprintf(fid, '%s\n', mkstr(context, ...
    '/! Z  CP  CS  RHO(g/cm3)  ATTN_P  ATTN_S', 57));

%% Others

% % min/max phase speed
% fprintf(fid, mkstr(sprintf('%.3f  %.3f', [PHASEMIN PHASEMAX]), '! min/max phase speed\n',57));

% % rmax
% fprintf(fid, mkstr(sprintf('%.3f', RMAX), '! RMAX (km)\n', 57));

% src
% assert(length(src_depth) == 2);
% fprintf(fid, '%d  /! NSD\n', nsd);
% fprintf(fid,'%.3f  %.3f / ! source depth \n', src_depth);
if length(src_depth)==1
    fprintf(fid, '%d  /! NSD\n', nsd);
    fprintf(fid,'%.3f  / ! source depth \n', src_depth);
else
    assert(length(src_depth) == 2 );
    fprintf(fid, '%d  /! NSD\n', nsd);
    fprintf(fid,'%.3f  %.3f / ! source depth \n', src_depth);
    
end

% rcv
assert(length(rcv_depth) == 2);
fprintf(fid,'%d  /! NRD \n', nrd);
fprintf(fid,'%.3f  %.3f / ! receiver depth \n',rcv_depth);

% NR
fprintf(fid,'%d  /! NR\n', length(Pos.r.range));
fprintf(fid,'%d %d  /! R(1:NR)\n', 0,RMAX);
% Option
context = sprintf('''%s''', Option);
fprintf(fid, '%s\n', mkstr(context, '/! R/C/I/S', 57));

% fprintf(fid,'%s\n', Option);
% Nbeams
fprintf(fid,'%d\n', Nbeams);
% ALPHA
assert(length(alpha) == 2);
fprintf(fid,'%.3f  %.3f / ! ALPHA1, 2 (degrees) \n', alpha);
% step
fprintf(fid,'%d %d %d / ! STEP (m), ZBOX (m), RBOX (km)',0,MEDIA{1}.zmax,RMAX);


fclose(fid);

end

function str = mkstr(str1, str2, str2_start)
len = length(str1);
while str2_start <= len
    str2_start = str2_start + 4;
end
dummy = char(32*ones(1,str2_start-len-1));
str = [str1, dummy, str2];
end