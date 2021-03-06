% MATLAB Code

% Var DATAs
b = [ 0 1 0 ] ;
pns = randi ( [ 0 1 ] , 1 , 20 ) ;
sn = length ( b ) ; 
pb = length ( pns ) ;
tb = sn * pb;
d = kron ( b , ones ( 1 , length ( pns ) ) ) ;
pn = reshape ( kron ( pns' , ones( 1 , sn ) ) , 1 , tb ) ;
tx = xor ( d , pn ) ;
bp = 2 * d - 1 ;
bpd = 2 * tx - 1 ;

% Baseband Signal
T = 1e-5 ;
Ns = 50 * pb ;
dt = T / Ns ;
t = dt : dt : T * sn ;
a = kron ( bp , ones ( 1 , 50 ) ) ;
ad = kron ( bpd , ones ( 1 , 50 ) ) ;

% Noise Channal
n = 2 * ( rand ( 1 , length ( ad ) ) > 0.25 ) - 1 ;
ad = ad .* n ;

% BPSK Modulation
fs = 1 / dt ;
fc = 1e6 ; wc = 2 * pi * fc ; 
s = a .* cos ( wc * t ) ;
sd = ad .* cos ( wc * t ) ;

% Orthogonal condition - signal detection
r = sum ( reshape ( sd .* cos ( wc * t ) , 50 , tb ) ) ;

for m = 1 : tb
    if r ( m ) >= 0 r ( m ) = 1 ;
    else r ( m ) = 0 ; end
end
    
rx = xor ( r , pn );
base = kron ( d , ones ( 1 , 50 ) );
base_pn = kron ( pn , ones ( 1 , 50 ) );
base_tx = kron ( tx , ones ( 1 , 50 ) );
base_rx = kron ( rx , ones ( 1 , 50 ) );

disp( '>> SND Data' ) , b
disp( '>> Symb Length PNCode' ) , pns
disp( '>> RECV Code Sample ReArrange' ) , rxs = reshape ( rx' , pb , sn )
disp( '>> RECV Restored Data' ) , d = unique ( rxs , 'rows' )

% Signal
BW = fs / 2 ;
f = linspace ( -BW , BW , length ( t ) ) ;
F = fftshift ( fft ( s ) ) / fs ;
Fd = fftshift ( fft ( sd ) ) / fs ;

% DRCT DIFS - CDMA
figure ( 'DRCT DIFS - CDMA' )
subplot ( 8 , 1 , 1 , 'align' )
plot ( t , base , 'color' , [ 0 0.7 0.4 ] , 'linewidth' , 2 )
axis ( [ min ( t ) max ( t ) -0.5 1.5 ] ) , ylabel ( 'DATA' )

title ( 'DS-BPSK' )
subplot ( 8 , 1 , 2 , 'align' )
plot ( t , base_pn , 'color' , [ 0 0.6 1 ] , 'linewidth' , 2 )
axis ( [ min ( t ) max ( t ) -0.5 1.5 ] ) , ylabel ( 'PN' )

subplot ( 8 , 1 , 3 , 'align' )
plot ( t , base_tx , 'm' , 'linewidth' , 2 )
axis ( [ min ( t ) max ( t ) -0.5 1.5 ] ) , ylabel ( 'Tx' )

subplot ( 8 , 1 , 4 , 'align' )
plot ( t , base_rx , 'g' , 'linewidth' , 2 )
axis ( [ min ( t ) max ( t ) -0.5 1.5 ] ) , ylabel ( 'Rx' )

i = find ( abs ( f ) <= 0.2 * BW ) ;
subplot ( 4 , 1 , 3 , 'align' )
mF = abs ( F ) ; plot ( f ( i )  , mF ( i ) , 'color' , [ 0 0.7 0.4 ] )
axis ( [ min ( f ( i ) ) max ( f ( i ) ) min ( mF ( i ) ) max ( mF ( i ) )
ylabel ( '|S(f)|' ) , grid
subplot( 4 , 1 , 4 , 'align' )
mFd = abs ( Fd ) ;
plot ( f ( i ) , mFd ( i ) , 'r' )
axis ( [ min ( f ( i ) ) max ( f ( i ) ) min ( mF ( i ) ) max ( mF ( i ) ) ] )
xlabel ( 'f [ HZ ]' ) , ylabel ( '|S_{DS}(F)|' ) , grid
