function [itr,src_func] = preprocrf(rf0,param)

TIME = param.Ti;
t1 = -3;
t2 = 20;
dt = param.dt;
nt1 = abs(t1)/dt;
nt2 = t2/dt;
x = size(rf0,2);

% normaliztion
itr = rf0./max(rf0(:));

% extract wavelet
srctmp=mean(itr,2);
[win] = waveform_win(srctmp,TIME,t1,1,-1);
src=srctmp.*win;
src_func = src./max(src(:));

% taper RF to remove later conversions
[win] = waveform_win(srctmp,TIME,t1,t2-5,3);
win = win*ones(1,size(itr,2));
itr = itr.*win;

% src_func = gradient(src_func,-0.5);
figure
wigb(itr(nt1:nt2,:),1.2,1:x,TIME(nt1:nt2))

end