function [plv, f] = cplxplv(x,params)
% Function to calculate plv using complex PCA and multitaper method
%
% USAGE: [plv,f] = cplxplv(x,params);
%
% x: input in channels x time x trials or time x trials
% params: A structure as in mtspectrumc() of the chronux toolbox i.e
%  params.Fs - Sampling rate
%  params.fpass - Frequency range of interest
%  params.pad - To zero pad or not to the next power of 2
%  params.tapers - [TW K]  i.e [time-half-BW product, #tapers]
%
%  plv: mean plv across tapers
if nargin<2
    params.Fs = 4096;
    params.fpass = [0 1000];
    params.pad = 0;
    params.tapers = [1 1];
end

if (ndims(x) == 3)
    nchans = size(x,1);
    ntime = size(x,2);
    ntrials = size(x,3);
else
    if(ndims(x) == 2)
        nchans = 1;
        ntime = size(x,1);
        ntrials = size(x,2);
        x = reshape(x,[1 ntime ntrials]);
    else
        fprintf(2,'\nSorry! The input data should be 2 or 3 dimensional!\n');
        help cplxplv
        return;
    end
end

ntaps = params.tapers(2);
TW = params.tapers(1);
w = dpss(ntime,TW,ntaps);

if(params.pad == 1)
    nfft = 2^nextpow2(ntime);
else
    nfft = ntime;
end

plv = zeros(ntaps,nfft);

x_dc = mean(x,2);
x = x - repmat(x_dc,[1,size(x,2),1]);

for k = 1:ntaps
    wbig = repmat(reshape(w(:,k),1,numel(w(:,k))),[nchans 1 ntrials]);
    X = fft(wbig.*x,nfft,2);

    Nf = nfft;
    spec = zeros(Nf,1);

    Csd = zeros(nchans,nchans,Nf); % Cross-spectral density matrix
    V = zeros(nchans,nchans,Nf); % Eigen vectors of the CSD matrix
    wtsc = zeros(nchans,Nf);
    
    X_norm = squeeze(mean(X./abs(X),3));
    
    for j = 1:Nf;
%         fprintf(1,'Doing frequency %d of %d\n',j,Nf)
%         Csd(:,:,j) = X_norm(:,j)*X_norm(:,j)'; % *
%         Csd(:,:,j) = (tempCsd + tempCsd')/2;
%         [V(:,:,j), D] = eig(squeeze(Csd(:,:,j))); % *
        for m = 1:nchans
            for n = 1:nchans
                CSD(m,n) = X_norm(m,j).*conj(X_norm(n,j));
            end
        end
        D = eig(CSD);
%         D = eig(X_norm(:,j)*X_norm(:,j)');
%         wtsc(:,j) = V(:,end,j); % Leading eigen vector *
        spec(j) = abs(real(D(end,end))/nchans);
        spec(j) = D(end,end)/nchans;
        %    wtsc(:,j) = wtsc(:,j).*spec(:,j);
    end
    plv(k,:) = spec;    
end

plv = squeeze(mean(plv,1));
f = (0:(nfft-1))*params.Fs/nfft;

plv = plv(:,(f >= params.fpass(1))&(f <= params.fpass(2)));
f = f((f >= params.fpass(1))&(f <= params.fpass(2)));

    
    
    
