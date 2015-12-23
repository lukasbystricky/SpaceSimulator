function srt = calculate_sidereal_time(date)

%omega = 7.292115856e-5;
omega = 1.0027379093*2*pi;
theta0 = 1.74933340;
% 
%t = 86636.55536351999*(datenum(date) - datenum('01-01-1970'));
t = (datenum(date) - datenum('01-01-1970'));

srt = wrapTo2Pi(t*omega + theta0);