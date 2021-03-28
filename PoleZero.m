clc;close all;
t=0.25;
fo=1/t;
Wo=2*pi*fo;
B=1.5; %damping
A=2080;

%T(s)=(�(s))/(G(s))=s^2/(s^2+2BWos+Wo^2 )
polinom1=[1 0 0];
polinom2=[1 2*B*Wo Wo^2];
Z=roots(polinom1);
P=roots(polinom2);


[numa,dena]=zp2tf(Z,P,A); %polezero to TF
F1=10^-2;
F2=10^2;

lg=logspace(-2, 2, 8192); %freqs range
lg=lg*(2*pi);

[R,F]=freqs(numa,dena,lg);  		%calculate frequency response
Mutlak=abs(R);
lh=semilogx((F./(2*pi)),Mutlak,'r');  %plot it (for impulse response)
figure;
lh=semilogx(F./(2*pi),20*log10(Mutlak/max(Mutlak)),'r');  %plot it for frequency response

xlabel('Frequency (Hz)')
ylabel('Magnitude (dB)')
title('Frequencey Response plot for:')

x=[F1 F2];
y=[-3 -3];

lh=line(x,y);
set(lh,'LineStyle',':');


% % % % a=1; 
% % % % h=0.707;        %damping
% % % % T=1;          %period
% % % % W0=2*pi*(1/T);  %omega

% % % % b=2*h*W0;       
% % % % c=W0^2;

% % % % polinom=[a b c];        % as^2+bs+c
% % % % sonuc=roots(polinom);   %radian
% % % % % Solve the equation  $x^4 - 1 = 0$.

% % % % % Create a vector to represent the polynomial, then find the roots.

% % % % % % % % p = [1 0 0 0 -1];
% % % % % % % % r = roots(p)

% % % % sonuc1=sonuc./(2*pi)    %frekans


% % % % a=1;
% % % % h=0.707;
% % % % f=50;
% % % % W0=2*pi*f;

% % % % b=2*h*W0;
% % % % c=W0^2;

% % % % polinom=[a b c];
% % % % sonuc=roots(polinom);
% % % % sonuc2=sonuc./(2*pi)

% % % %  Z=[0];
% % % %         P=[sonuc1; sonuc2];
% % % %         A=[50]

% % % % %if isempty(Z),
% % % % %    error(['Type "',TYPE,'" is not a recognised sensor type']);
% % % % %    return;
% % % % %end

% % % % Z=((Z)*(2*pi));             %convert to Radians / sec
% % % % P=((P)*(2*pi));             %convert to Radians / sec
% % % % A=A* ((((2*pi)^(length(P))) / ((2*pi)^(length(Z)))));

% % % % z=P;
% % % % p=Z;
% % % % k=A;
% % % % plotpz_blue(p,z,k)
