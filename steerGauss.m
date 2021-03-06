function R = steerGauss(I,Wsize,sigma,theta) 
% 将角度转化在[0,pi]之间
theta = theta/180*pi;
% 计算二维高斯核在x,y方向的偏导gx,gy
k = [-Wsize:Wsize];
g = exp(-(k.^2)/(2*sigma^2));
gp = -(k/sigma).*exp(-(k.^2)/(2*sigma^2));
gx = g'*gp;
gy = gp'*g;
% 计算图像Ｉ在x,y方向的滤波结果
Ix = conv2(I,gx,'same');
Iy = conv2(I,gy,'same');
% 计算图像Ｉ在theta方向的滤波结果 axis image; colormap(gray);
J = cos(theta)*Ix+sin(theta)*Iy;
figure,imshow(J);
figure,subplot(1,3,1),imshow(I),title('原图像');
subplot(1,3,2),axis image;imshow(cos(theta)*gx+sin(theta)*gy),title('滤波模板');
subplot(1,3,3),axis image;imshow(J),title('滤波结果'); 