
Nimg = 6;
fbase = '../../../../Desktop/umpc_mm_data/trial36_';

Aout = imgin([fbase,int2str(1),'.png'], Nimg);

for i=2:Nimg
	A = imgin([fbase,int2str(i),'.png'], Nimg);
	Aout = imfuse(A, Aout, 'method', 'blend','Scaling','joint');
% 	Aout = A+Aout;
end

% 	Aout = imadjust(Aout);
imshow(Aout)
imcontrast

function B=imgin(fname, Nimg)
	B = imread(fname);
	B = rgb2gray(B);
	B = imadjust(B)/Nimg;
end
