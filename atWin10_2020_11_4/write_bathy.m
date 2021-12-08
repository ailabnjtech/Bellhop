function write_bathy( bathyfil, Option, bathy )

% Write a bathy file

if ( ~strcmp( bathyfil( end - 3 : end ), '.bty' ) )
  bathyfil = [ bathyfil '.bty' ]; % append extension
end

fid = fopen( bathyfil, 'w' );

fprintf(fid, '%s\n', mkstr(sprintf('''%s''', Option), '! Option', 65));

fprintf(fid, '%d\n', length(bathy));

   for ii = 1:length(bathy)
        fprintf(fid, '%8.4f %8.4f /\n', bathy(ii,1),bathy(ii,2));
    end






fclose( fid );

end

function str = mkstr(str1, str2, str2_start)
len = length(str1);
while str2_start <= len
    str2_start = str2_start + 4;
end
dummy = char(32*ones(1,str2_start-len-1));
str = [str1, dummy, str2];
end