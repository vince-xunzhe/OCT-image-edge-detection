clear all;
close all;
clc;

img=imread('OCT1.jpg');
imshow(img);title('');
%img=double(img);

img=rgb2gray(img);
%img=double(gray);
[m n]=size(img);
%img_canny=edge(g,'canny',[0.03,0.06]); 
%subplot(222);imshow(img_canny);title('Canny');
%st=strel('square',10);
%img=imclose(o_img,st);
%img=rgb2gray(img);
%img=im2double(img);

img=wiener2(img,[10 10]);
f=img;
%img=filter2(fspecial('average',10),img);
%img=medfilt2(img,[8 8]);

figure,imshow(img);title('');

%%高斯滤波
w=fspecial('gaussian',[5 5]);    
img=imfilter(img,w,'replicate');
figure,imshow(uint8(img));title('高斯滤波');

%%sobel边缘检测
w=fspecial('sobel');
img=double(img);
img_w=imfilter(img,w,'replicate');                                         %求横边缘
w=w';
img_h=imfilter(img,w,'replicate');                                         %求竖边缘
img=sqrt(img_w.^2+img_h.^2);                                               %平方和在开方
figure,imshow(uint8(img));

%%下面是非极大抑制
new_edge=zeros(m,n);
for i=2:m-1
    for j=2:n-1
        Mx=img_w(i,j);
        My=img_h(i,j);
        
        if My~=0
            o=atan(Mx/My);      %边缘的法线弧度
        elseif My==0 && Mx>0
            o=pi/2;
        else
            o=-pi/2;            
        end
        
        %Mx处用My和img进行插值
        adds=get_coords(o);                                                %边缘像素法线一侧求得的两点坐标，插值需要       
        M1=My*img(i+adds(2),j+adds(1))+(Mx-My)*img(i+adds(4),j+adds(3));   %插值后得到的像素，用此像素和当前像素比较 
        adds=get_coords(o+pi);                                             %边缘法线另一侧求得的两点坐标，插值需要
        M2=My*img(i+adds(2),j+adds(1))+(Mx-My)*img(i+adds(4),j+adds(3));   %另一侧插值得到的像素，同样和当前像素比较
        
        isbigger=(Mx*img(i,j)>M1)*(Mx*img(i,j)>=M2)+(Mx*img(i,j)<M1)*(Mx*img(i,j)<=M2); %如果当前点比两边点都大置1
        
        if isbigger
           new_edge(i,j)=img(i,j); 
        end        
    end
end
figure,imshow(uint8(new_edge));


%%下面是滞后阈值处理
up=100;                                                                    %上阈值 216 140 25
low=70;                                                                    %下阈值
set(0,'RecursionLimit',10000);                                             %设置最大递归深度
for i=1:m
    for j=1:n
      if new_edge(i,j)>up &&new_edge(i,j)~=255                             %判断上阈值
            new_edge(i,j)=255;
            new_edge=connect(new_edge,i,j,low);
      end
    end
end
figure,imshow(new_edge==255);

Dx=[1,0,-1];
Dy=[1;0;-1];
b=imfilter(f,Dx);
b=imfilter(b,Dy);
%b=gradient(double(f));
figure,imshow(b);title('轴向梯度信息');

vector_1=zeros(m);
vector_2=zeros(m);
vector_1=f(m/2,1:n);
f2=uint8(new_edge);
vector_2=b(m/2,1:n);
figure,stem(vector_1);
figure,stem(vector_2);


%x=b.*(uint8(new_edge));
x=double(b).*(new_edge==255);
figure,imshow(x);title('梯度信息与边缘结果');
%plate=zeros(m,n);
%plate_1=zeros(m,n);
%for i=1:m
%    for j=1:n
%        plate(i,j)=new_edge(i,j);
%        plate_1(i,j)=new_edge(i,j);
%    end
%end
%plate=double(plate);
%plate_1=plate.*f;
%figure,subplot(1,2,1);imshow(plate);
%subplot(1,2,2);imshow(plate_1);

%steerGauss(uint8(new_edge),7,1,45);


sigma=0.5;
h=floor(2*sigma+1); 
w=h; 
[x y]=meshgrid(-w:w,-h:h);
Ga=0.9213*(-2.254*x+x.^3).*exp(-(x.^2+y.^2)/(2*sigma^2));   %各种三阶的幅度系数
Gb=1.843*(-0.7515+x.^2).*y.*exp(-(x.^2+y.^2)/(2*sigma^2));
Gc=0.9780*(-0.7515+y.^2).*x.*exp(-(x.^2+y.^2)/(2*sigma^2));
Gd=0.9780*(-2.254*y+y.^3).*exp(-(x.^2+y.^2)/(2*sigma^2));

img=double(f);
[m n]=size(img);
edge=zeros(m,n);

for i=0:30:360              %一次转过30度角
    theta=(i/180)*pi;
    Ka=cos(theta)^3;            %各种三阶的角度系数
    Kb=-3*cos(theta)^2*sin(theta);
    Kc=3*cos(theta)*sin(theta)^2;
    Kd=-sin(theta)^3;
    G=Ka*Ga+Kb*Gb+Kc*Gc+Kd*Gd;      %待卷积模板
    
    imgn=imfilter(img,G,'replicate');
    
    figure,imshow(imgn,[]);
    edge=sqrt(edge.^2+imgn.^2);
end

figure,imshow(edge,[])











    




