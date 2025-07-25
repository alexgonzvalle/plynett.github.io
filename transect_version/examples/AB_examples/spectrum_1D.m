
function [amps,periods,phases] = spectrum_1D(f,h,Hs,per,H_IG,T_IG);

f_peak=1/per;    % Frequency here is defined as 1/T - this is NOT angular frequency (Hz)
nf = length(f);

% Shallow water TMA spectrum
g=9.81;		% gravity
gamma=3.0;
beta=0.0624/(0.23+0.033*gamma-0.185*(1.9+gamma)^-1);

for i=1:length(f)
    omega_h=2*3.1415*f(i)*sqrt(h/g);
    if omega_h>2 % Note that adding this phiK term to the JONSWAP spectrum as given here will
        % cause the intergrated E, as found below, to no longer equal Hs.  This will need to be corrected
        phiK=1;
    elseif omega_h<1
        phiK=0.5*omega_h^2;
    else
        phiK=1-0.5*(2-omega_h)^2;
    end

    if f(i)<=f_peak
        sigma=0.07;
    else
        sigma=0.09;
    end

    frat=f(i)/f_peak;

    E(i)=beta*Hs^2/(f(i)*frat^4)*exp(-1.25/frat^4)*gamma^(exp(-(frat-1)^2/(2*sigma^2)))*phiK;  % energy density
end


Hmo = 0;

for i=1:length(f)
    if i==1
        del_f=(f(2)-f(1));
    elseif i==length(f)
        del_f=(f(length(f))-f(length(f)-1));
    else
        del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
    end

    Hmo = Hmo + E(i)*del_f;
end

Hmo_full_spectrum = sqrt(Hmo)*4.004;

% make phiK correction to keep Hs correct
if Hmo_full_spectrum>0
    E=E*(Hs/Hmo_full_spectrum)^2;
end


% Truncate ends of spectrum at values of 5% of the max,

Hmo_truncated_spectrum=0;
amps=zeros(nf,1);
periods=amps;
phases=amps;
for i=1:nf
    if i==1
        del_f=(f(2)-f(1));
    elseif i==nf
        del_f=(f(nf)-f(nf-1));
    else
        del_f=(f(i+1)-f(i))/2.+(f(i)-f(i-1))/2.;
    end

    Hmo_truncated_spectrum = Hmo_truncated_spectrum + E(i)*del_f;
    amps(i)=sqrt(2*E(i)*del_f);
    periods(i) = 1/f(i);
    phases(i) = rand*2*3.1415;
end

amps = [H_IG/2; amps];
periods = [T_IG; periods];
phases = [0; phases];






