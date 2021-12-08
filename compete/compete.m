load('parameter20.mat')
load('data1_9181_9271.mat')

tl_interp = interp1(rv,20*log10(intMid),range);

cp = mse(tl_interp-(p_mf-200));

fprintf('cp=%f',cp)
