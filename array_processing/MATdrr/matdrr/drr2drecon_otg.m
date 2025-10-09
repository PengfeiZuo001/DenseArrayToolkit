function [ D1, D2 ] = drr2drecon_otg(D,x,nx,ox,mx,flow,fhigh,dt,N,K,Niter,eps,verb,mode,a)
%  DRR2DRECON_OTG: OTG DRR for 2D seismic reconstruction 
%
%  IN   D:   	 intput 2D data [time x traces]
%       x:     input x coordinates [1 x ntraces]
%       nx:    input number of binned x points
%       ox:    min of x
%       mx:    max of x
%       flow:   processing frequency range (lower)
%       fhigh:  processing frequency range (higher)
%       dt:     temporal sampling interval
%       N:      number of singular value to be preserved
%       K:     damping factor (default: 4)
%       Niter:  number of maximum iteration
%       eps:    tolerence (||S(n)-S(n-1)||_F<eps)
%       verb:   verbosity flag (default: 0)
%       mode:   mode=1: denoising and reconstruction
%               mode=0: reconstruction only
%       a:      weight vector
%
%  OUT  D1:  	output data on regular grid [time x nx]
%       D2:     output data on original grid [time x ntraces]
%
%  Copyright (C) 2021 The University of Texas at Austin
%  Copyright (C) 2021 Yangkang Chen
%  Adapted for 2D by Yunfeng Chen, 2025
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%
%  References:
%
%  [1] Chen, Y., W. Huang, D. Zhang, W. Chen, 2016, An open-source matlab code package for improved rank-reduction 3D seismic data denoising and reconstruction, Computers & Geosciences, 95, 59-66.
%  [2] Chen, Y., D. Zhang, Z. Jin, X. Chen, S. Zu, W. Huang, and S. Gan, 2016, Simultaneous denoising and reconstruction of 5D seismic data via damped rank-reduction method, Geophysical Journal International, 206, 1695-1717.
%  [3] Huang, W., R. Wang, Y. Chen, H. Li, and S. Gan, 2016, Damped multichannel singular spectrum analysis for 3D random noise attenuation, Geophysics, 81, V261-V270.

if nargin==0
    error('Input data must be provided!');
end

if nargin==1
    error('Sampling coordinates should be given');
end

if nargin==2
    flow=1;
    fhigh=124;
    dt=0.004;
    N=1;
    K=4;
    Niter=30;
    eps=0.00001;
    verb=0;
    mode=1;
end;

if mode==0;
    a=ones(1,Niter);
end

nt=size(D,1);
D1=zeros(nt,nx);

nf=2^nextpow2(nt);

% Transform into F-X domain
DATA_FX=fft(D,nf,1);
DATA_FX0=zeros(nf,nx);

% First and last nts of the DFT.
ilow  = floor(flow*dt*nf)+1;

if ilow<1;
    ilow=1;
end;

ihigh = floor(fhigh*dt*nf)+1;

if ihigh > floor(nf/2)+1;
    ihigh=floor(nf/2)+1;
end

lx=floor(nx/2)+1;
lxx=nx-lx+1;
M=zeros(lx,lxx);

% construct par
par.x=x;
par.nx=nx;
par.ox=ox;
par.mx=mx;
s=1*ones(Niter,1);

% main loop
for k=ilow:ihigh
    S_obs=squeeze(DATA_FX(k,:)).'; %1D vector  
    Sn_1=zeros(nx,1);
    
    for iter=1:Niter
        Sn=Sn_1-s(iter)*inter_op_1d(inter_op_1d(Sn_1,par,-1)-S_obs,par,1);
        
        M=P_H_1d(Sn,lx);
        M=P_RD(M,N,K);
        Sn=P_A_1d(M,nx,lx);
        
        if norm(Sn-Sn_1,'fro')<eps
            break;
        end
        Sn_1=Sn;
    end
    
    DATA_FX0(k,:) = DATA_FX0(k,:)+reshape(Sn,1,nx);
    
    if(mod(k,5)==0 && verb==1)
        fprintf( 'F %d is done!\n\n',k);
    end
    
end

% Honor symmetries
for k=nf/2+2:nf
    DATA_FX0(k,:) = conj(DATA_FX0(nf-k+2,:));
end

% Back to TX (the output)
D1=real(ifft(DATA_FX0,[],1));
D1=D1(1:nt,:);
D2 = zeros(size(D));
for n = 1:size(D1,1)
    D2(n,:)=inter_op_1d(squeeze(D1(n,:)),par,-1);
end
return

function [dout]=P_H_1d(din,lx)
% forming Hankel matrix for 1D
% din: [nx x 1]
nx=length(din);
lxx=nx-lx+1;
dout=hankel(din(1:lx),[din(lx:nx)]);
return

function [dout]=P_RD(din,N,K)
% Rank reduction on the Hankel matrix
[U, S, V] = svds(din, N+1);  % Compute first N+1 singular values
S = diag(S);
S(1:N) = S(1:N) .* (1 - S(N+1)^K ./ (S(1:N).^K + eps));
dout = U(:,1:N) * diag(S(1:N)) * V(:,1:N)';
return

function [dout]=P_A_1d(din,nx,lx)
% Averaging the Hankel matrix to output the result
lxx=nx-lx+1;
dout=zeros(nx,1);

for i=1:nx
    if i<lx
        count = i;
        for j=1:i
            dout(i) = dout(i) + din(j, i-j+1) / count;
        end
    else
        count = nx - i + 1;
        for j=1:count
            dout(i) = dout(i) + din(lx-j+1, i-lx+j) / count;
        end
    end
end
return

function [U] = inter_op_1d(D,par,adj)
% 1D interpolation operator
% if adj = 1, D: 1D data (original traces), U: 1D data (gridded)
% if adj = -1, D: 1D data (gridded), U: 1D data (original traces)

rx=par.x;
dx=(par.mx-par.ox)/(par.nx-1);
xx=par.ox+[0:par.nx-1]*dx;

Nu=length(rx);
Nx=length(xx);

if adj==1
    U=zeros(Nx,1);
else
    U=zeros(Nu,1);
end

for k=1:Nu
    ia=floor((rx(k)-xx(1))/dx)+1;
    ib=ia+1;
    
    if ib>1 && ib<=Nx
        t=(rx(k)-xx(ia))/dx;
        if adj==1
            U(ia)=U(ia)+(1-t)*D(k);
            U(ib)=U(ib)+t*D(k);
        else
            U(k)=(1-t)*D(ia)+t*D(ib);
        end
    end
end

return
