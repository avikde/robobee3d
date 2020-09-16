function wltest()

clear mex

% repeatedly call for a made-up pdes to test
Nt = 500;
U = zeros(Nt,4);
pdotdes = [0,0,10,0,0,0];

U(1,:) = [140.0, 0., 0., 0.];
mb = 100;
g = 9.81e-3;
h0 = [0, 0, mb * g, 0, 0, 0];

% poptsEmp2
popts = [1.1194660044615514 , -0.012509989456907032 , -1.0861614307812186 , 2.4851677856988856 , 1.0 , 6.844974694727402e-05 , 0.006913760202682983 , -0.014818711614996979 , 1.0 , 0.7099927690180734 , -0.06888578769654327 , 1.0 , 2.2007661373369336 , 
1.0 , 1.0 , -0.0020225622871363996 , 1.0273004470130044e-05 , 0.09269580456410237 , 0.036251736135975196 , 1.0 , 8.734129661583814e-09 , -0.0006201755761980693 , -0.00017875339271336048 , 1.0 , -0.05034637957148408 , -0.040900504208754775 , 1.0 , -0.0007749399045436298 , 1.0 , 1.0 , -1.210928046887087 , 0.003924244409026407 , 1.1473139983862783 , -10.093444365164013 , 1.0 , 0.0001219293546880897 , -0.01106560071988489 , 0.04971229082512828 , 1.0 , -10.182796228428456 , 0.7704674202119165 , 1.0 , 3.7037292384948968 
, 1.0 , 1.0 , -58.54651372841956 , 0.6426626139787593 , 5.896120168505504 , 27.068153274917073 , 1.0 , -0.0032495300808689518 , -0.055034724705108805 , -0.40638954827328905 , 1.0 , 8.366443418119623 , 8.999567348553631 , 1.0 , -104.22344178880671 , 1.0 , 1.0 , -6.7303575189310765 , 0.0381719425810673 , 17.256396078192687 , -2.3122815061661517 , 1.0 , -0.00018761860271511883 , -0.12584828582103202 , 0.01070424287314799 , 1.0 , 9.982972380159119 , 2.1298243744292846 , 1.0 , -5.831985522884596 , 1.0 , 1.0 , -12.657677705070068 , 0.14703931014183502 , 0.9962740096920417 , -24.193455539723775 , 1.0 , -0.0008767858139277795 , -0.011349289130157127 
, 0.13217466383249077 , 1.0 , 9.75767475739342 , 12.111853855954712 , 1.0 , -27.971730203492886 , 1.0 , 1.0];
% % poptsEmp row major
% popts = [-0.41498716334913877, 0.004161530470390393, -0.04921896208897451, -0.08114084916912893, 1.0, -2.1498519749986776e-05, 0.0007452193674740687, -0.000592201745189127, 1.0, 2.489837075232938, -0.8220389384978928, 1.0, -0.6424779320800695, 1.0, 1.0, 0.04853925386307669, -0.0005410195533586718, 0.03325717298367986, 0.11402108349065512, 1.0, 2.9591070779816484e-06, -0.0003656669558602931, -0.0006235124178038529, 1.0, -0.15022060357201647, 0.0003189736297552722, 1.0, 0.09522768136462195, 1.0, 1.0, -1.0443154298902466, 0.008730954429386921, 3.774256296882532, -6.565041851772019, 1.0, 7.255944919680254e-05, -0.024440375440139876, 0.03221142824195457, 1.0, -23.321197610424026, 3.865481128472532, 1.0, 1.7626691411728337, 1.0, 1.0, -39.162566526689496, 0.4347272529556003, 35.20801123124533, 8.960105550238962, 1.0, -0.002179040590683337, -0.20825719312770322, -0.2710167102094739, 1.0, 98.4856348629078, 34.454898506387536, 1.0, -68.47304863790173, 1.0, 1.0, 3.86697945239634, -0.054406534814573505, 4.136310537263415, -2.40430210592086, 1.0, 0.00034641731566390497, -0.05418716099028047, 0.011121246757241961, 1.0, 5.142757475687521, -10.274310705200136, 1.0, 11.418953851491452, 1.0, 1.0, 1.5868629373519454, -0.013584604586031352, 8.897152221477027, 6.856427471376807, 1.0, 1.2282136287999407e-05, -0.05776369789309159, -0.0360868798609794, 1.0, 8.468979086282319, 0.7365484310116243, 1.0, 0.44877875839694553, 1.0, 1.0];
% %popts (numerical from simulation)
% popts = [-0.006546634426549249, 5.619730433095801e-05, -0.10046126965362534, 2.8269636386847336e-07, 3.3260105206082036e-07, -8.561129105121841e-07, 0.0010378391603621869, -2.1526976705866592e-09, -2.4682847177637722e-09, 0.00030895094854226844, -5.498973958767363e-08, -3.6850719187462874e-07, -0.012527592056397634, -0.2386723461724283, -0.011289258444268205, 9.029009225475524e-06, -1.4986231356039162e-07, 4.371403330671029e-08, 1.023574109414231, 0.017039273759093976, 1.2043402924645282e-09, -3.473465666769341e-10, -0.012121931849592966, -0.00020828441245571597, 4.490548891202895e-09, -0.0028584701763253568, 0.00042179444241666956, -8.389602378310228e-06, -8.792969701568155e-07, 7.233850163004005e-06, -0.4235793432635114, 0.0021338360840689285, 2.4430560265465476e-05, 3.328050926551559e-05, 1.0836641499138236e-05, 0.00023127061314787283, -1.9501837214353561e-07, -2.513252554856101e-07, -7.917131239703786e-08, 1.911118527187613e-05, -2.529699081072906e-05, -3.1099855170076864e-05, 2.4948092705639127, -0.007818245388260668, 0.10682098780722615, 0.0004970273589210487, -6.882776500900782e-06, -0.0004592020239540682, -28.07432375743903, -0.0008354886614699838, 5.190747088376899e-08, 3.322441340705155e-06, 0.521662596083329, 3.863531490483966e-05, -0.0004502665930561039, 0.01632414175861271, 0.01085839020243772, 5.8402721126669827e-05, -0.0013350535953576411, -0.003619350491662568, -0.04967855748617539, 0.0008594741997761248, -5.141024035634106, -3.6947357033333146e-05, 0.0001902321337739334, -9.815504649772921e-06, 0.06222625213186152, 2.897975880464687e-07, -1.5177479898802682e-06, 0.0011965904064848794, -8.836115186850822e-05, -6.286132183531711e-05, -0.16089773884138625, -0.5157669811698968, -0.017986323101450558, 2.7799826925561083e-05, -4.33980866463365e-07, -3.7886225940335855e-07, -0.30424162339652516, -3.041423746090091, 3.4860207993083713e-09, 7.844144603352391e-09, 0.0033650097914274413, 0.022997043353002334, -2.9579795830697566e-06, -0.4118396574989255, -0.006259667773324984, -0.00020662428061023186, 7.841637217987254e-06, 3.213543642962808e-05];


for i = 2:Nt
	if i > 100
			pdotdes(3) = 20;
	end
	if i > 200
			pdotdes(3) = 30;
	end
	if i > 300
			pdotdes(3) = 40;
	end
	if i > 400
			pdotdes(3) = 50;
	end
	U(i,:) = wlControllerUpdate(single(U(i-1,:)), single(h0), single(pdotdes), single(popts), 1000);
end

% convert
Vmean = U(:,1);
uoffs = -U(:,2);
udiff = U(:,3);

vleft = Vmean .* (1 + udiff) .* (1 + uoffs);
vright = Vmean .* (1 - udiff) .* (1 + uoffs);

drv_pch = Vmean .* uoffs;
drv_amp = abs((vleft-vright))/2+min(vleft,vright);
drv_roll = (vleft-vright)/4; 
drv_bias = max(vleft,vright)+2*abs(drv_pch); % voltage

subplot(2,1,1)
hold all
plot(Vmean)
plot(drv_bias)
plot(drv_amp)
ylim([100,220])
legend('Vmean','drv_bias','drv_amp')

subplot(2,1,2)
hold all
plot(drv_roll)
plot(drv_pch)
ylim([-30,30])
legend('drv_roll','drv_pch')

