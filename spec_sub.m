function ssdenoised_signal = spec_sub(noisy_signal)

winsize = 256;% window length

n = estimatenoise(noisy_signal);% noise level
 
fs = 4000;

size = length (noisy_signal);% speech length
numofwin = floor (size / winsize);% Number of window

% Hamming window defined
ham=hamming(winsize)';
hamwin=ones(1,size);
improved=zeros(1,size);
 
%Noisy signal
y = noisy_signal;


% Noise processing
%noisy=n*randn(1,winsize);
noisy = wgn(1,size, 10*log10(n),'real');
N=fft(noisy);

npow=abs(N);

for q=1:2*numofwin-1
yframe = y(1+ (q-1) * winsize / 2: winsize + (q-1) * winsize / 2);
yframe=yframe';
% framing hamwin(1+(q-1)*winsize/2:winsize+(q-1)*winsize/2)=hamwin(1+(q-1)*winsize/2:winsize+(q-1)*winsize/2)+ham;%
 
% Plus noise signal FFT
y1=fft(yframe.*ham);
ypow = abs(y1);% signal plus noise amplitude
yangle = angle(y1);% Phase
 

% Calculated power spectral density
Py=ypow.^2;
Pn=npow.^2;
%Pyy=ypow.^a;
%Pnn=npow.^a;
 
%The basic spectral subtraction%
for i=1:winsize
if Py(i)-Pn(i)>0
Ps(i)=Py(i)-Pn(i);
else
Ps(i)=0;
end
end
s=sqrt(Ps).*exp(j*yangle);

%{
for i=1:winsize
if Pyy(i)-b*Pnn(i)>0
Pss(i)=Pyy(i)-b*Pnn(i);
else
Pss(i)=c*Pnn(i);
end
end
ss=Pss.^(1/a).*exp(j*yangle);
%} 

% De-noising voice IFFT
temp2=improved(1+(q-1)*winsize/2:winsize+(q-1)*winsize/2)+real(ifft(s));
improved(1+(q-1)*winsize/2:winsize+(q-1)*winsize/2)=temp2;
end



% Removal of the gain due to Hamming window
for i=1:size
if hamwin(i)==0
improved(i)=0;
else
improved(i)=improved(i)/hamwin(i);
end
end

improved = improved';

ssdenoised_signal = improved;
