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

%%��˹�˲�
w=fspecial('gaussian',[5 5]);    
img=imfilter(img,w,'replicate');
figure,imshow(uint8(img));title('��˹�˲�');

%%sobel��Ե���
w=fspecial('sobel');
img=double(img);
img_w=imfilter(img,w,'replicate');                                         %����Ե
w=w';
img_h=imfilter(img,w,'replicate');                                         %������Ե
img=sqrt(img_w.^2+img_h.^2);                                               %ƽ�����ڿ���
figure,imshow(uint8(img));

%%�����ǷǼ�������
new_edge=zeros(m,n);
for i=2:m-1
    for j=2:n-1
        Mx=img_w(i,j);
        My=img_h(i,j);
        
        if My~=0
            o=atan(Mx/My);      %��Ե�ķ��߻���
        elseif My==0 && Mx>0
            o=pi/2;
        else
            o=-pi/2;            
        end
        
        %Mx����My��img���в�ֵ
        adds=get_coords(o);                                                %��Ե���ط���һ����õ��������꣬��ֵ��Ҫ       
        M1=My*img(i+adds(2),j+adds(1))+(Mx-My)*img(i+adds(4),j+adds(3));   %��ֵ��õ������أ��ô����غ͵�ǰ���رȽ� 
        adds=get_coords(o+pi);                                             %��Ե������һ����õ��������꣬��ֵ��Ҫ
        M2=My*img(i+adds(2),j+adds(1))+(Mx-My)*img(i+adds(4),j+adds(3));   %��һ���ֵ�õ������أ�ͬ���͵�ǰ���رȽ�
        
        isbigger=(Mx*img(i,j)>M1)*(Mx*img(i,j)>=M2)+(Mx*img(i,j)<M1)*(Mx*img(i,j)<=M2); %�����ǰ������ߵ㶼����1
        
        if isbigger
           new_edge(i,j)=img(i,j); 
        end        
    end
end
figure,imshow(uint8(new_edge));


%%�������ͺ���ֵ����
up=100;                                                                    %����ֵ 216 140 25
low=70;                                                                    %����ֵ
set(0,'RecursionLimit',10000);                                             %�������ݹ����
for i=1:m
    for j=1:n
      if new_edge(i,j)>up &&new_edge(i,j)~=255                             %�ж�����ֵ
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
figure,imshow(b);title('�����ݶ���Ϣ');

vector_1=zeros(m);
vector_2=zeros(m);
vector_1=f(m/2,1:n);
f2=uint8(new_edge);
vector_2=b(m/2,1:n);
figure,stem(vector_1);
figure,stem(vector_2);


%x=b.*(uint8(new_edge));
x=double(b).*(new_edge==255);
figure,imshow(x);title('�ݶ���Ϣ���Ե���');
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
Ga=0.9213*(-2.254*x+x.^3).*exp(-(x.^2+y.^2)/(2*sigma^2));   %�������׵ķ���ϵ��
Gb=1.843*(-0.7515+x.^2).*y.*exp(-(x.^2+y.^2)/(2*sigma^2));
Gc=0.9780*(-0.7515+y.^2).*x.*exp(-(x.^2+y.^2)/(2*sigma^2));
Gd=0.9780*(-2.254*y+y.^3).*exp(-(x.^2+y.^2)/(2*sigma^2));

img=double(f);
[m n]=size(img);
edge=zeros(m,n);

for i=0:30:360              %һ��ת��30�Ƚ�
    theta=(i/180)*pi;
    Ka=cos(theta)^3;            %�������׵ĽǶ�ϵ��
    Kb=-3*cos(theta)^2*sin(theta);
    Kc=3*cos(theta)*sin(theta)^2;
    Kd=-sin(theta)^3;
    G=Ka*Ga+Kb*Gb+Kc*Gc+Kd*Gd;      %�����ģ��
    
    imgn=imfilter(img,G,'replicate');
    
    figure,imshow(imgn,[]);
    edge=sqrt(edge.^2+imgn.^2);
end

figure,imshow(edge,[])











    




