function re=get_coords(angle)       %angle�Ǳ�Ե���߽Ƕȣ����ط���ǰ������
    sigma=0.000000001;
    x1=ceil(cos(angle+pi/8)*sqrt(2)-0.5-sigma);
    y1=ceil(-sin(angle-pi/8)*sqrt(2)-0.5-sigma);
    x2=ceil(cos(angle-pi/8)*sqrt(2)-0.5-sigma);
    y2=ceil(-sin(angle-pi/8)*sqrt(2)-0.5-sigma);
    re=[x1 y1 x2 y2];

end