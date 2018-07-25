function R = steerGauss(I,Wsize,sigma,theta) 
% ���Ƕ�ת����[0,pi]֮��
theta = theta/180*pi;
% �����ά��˹����x,y�����ƫ��gx,gy
k = [-Wsize:Wsize];
g = exp(-(k.^2)/(2*sigma^2));
gp = -(k/sigma).*exp(-(k.^2)/(2*sigma^2));
gx = g'*gp;
gy = gp'*g;
% ����ͼ�����x,y������˲����
Ix = conv2(I,gx,'same');
Iy = conv2(I,gy,'same');
% ����ͼ�����theta������˲���� axis image; colormap(gray);
J = cos(theta)*Ix+sin(theta)*Iy;
figure,imshow(J);
figure,subplot(1,3,1),imshow(I),title('ԭͼ��');
subplot(1,3,2),axis image;imshow(cos(theta)*gx+sin(theta)*gy),title('�˲�ģ��');
subplot(1,3,3),axis image;imshow(J),title('�˲����'); 