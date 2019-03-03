%初始化参数

a=3.174;
d=1.469;
e1=0;
e2=3;
e3=0;
e4=-3;
e5=-3;
e6=-3;
%以下参数中，下标a表示西格玛键，下标b表示派键，c表示德尔塔键
Vdda=-0.5;%调这个参数可以改变中间峰的高低，和两边的平缓程度
Vddb=1.5;%调这个参数可以改变边缘那一段的平缓和陡峭程度
Vddc=1.5;%调这个参数可以让K点附件的谷变深
Vpda=1.2;%调这个参数可以让e1那条带中间的峰变高
Vpdb=0.6;

Vppa=0;
Vppb=0;

%磁效应和SOC
c=0.3;
m1=0.3;
m2=0.3;

%画能带时的纵坐标、横坐标以及点密度
Band=[];
x_cordi=[];
accu=100;

%从M点到K+点
for order=0:accu
 %定义从M点到K+点的K值
  kx=-sqrt(3)/4-1/(4*sqrt(3))*order/accu;
  ky=1/4*(1-order/accu);
 %得到矩阵
  H=TBMatrix(a,d,e1,e2,e3,e4,e5,e6,Vdda,Vddb,Vddc,Vpda,Vpdb,Vppa,Vppb,c,m1,m2,kx,ky);
 %得到能带
  Band=[Band;eig(H)];
 %从M点到K+点的能带数据
  temp1=[order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2;order/2];
  x_cordi=[x_cordi;temp1];
end

%从K+点到K-点，途径Gamma点
for order=0:2*accu
 %从K+到K-点的K值
  ky=0;
  kx=-1/sqrt(3)+1/sqrt(3)*order/accu;
 %得到矩阵
  H=TBMatrix(a,d,e1,e2,e3,e4,e5,e6,Vdda,Vddb,Vddc,Vpda,Vpdb,Vppa,Vppb,c,m1,m2,kx,ky);
 %得到能带
  Band=[Band;eig(H)];
 %从K+到K-点的能带数据
  temp2=[accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order;accu/2+order];
  x_cordi=[x_cordi;temp2];
end

%从K-到M点
for order=0:accu
 %从K-到M点的K值
  kx=1/sqrt(3)-1/(4*sqrt(3))*order/accu;
  ky=1/4*order/accu;
 %得到矩阵
  H=TBMatrix(a,d,e1,e2,e3,e4,e5,e6,Vdda,Vddb,Vddc,Vpda,Vpdb,Vppa,Vppb,c,m1,m2,kx,ky);
 %得到能带
  Band=[Band;eig(H)];
 %从K-点到M点的能带数据
  temp3=[order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2;order/2+accu*5/2];
  x_cordi=[x_cordi;temp3];
end


%从M点到Gamma点的K值
%kx=sqrt(3)/4*(1-order/accu);
%ky=1/4*(1-order/accu);

%从K-点到M点的能带数据
%temp4=[order*sqrt(3)/2+accu*3,order*sqrt(3)/2+accu*3,order*sqrt(3)/2+accu*3,order*sqrt(3)/2+accu*3,order*sqrt(3)/2+accu*3,order*sqrt(3)/2+accu*3,order*sqrt(3)/2+accu*3,order*sqrt(3)/2+accu*3,order*sqrt(3)/2+accu*3];
%x_cordi=[x_cordi;temp4];

function [H]=TBMatrix(a,d,e1,e2,e3,e4,e5,e6,Vdda,Vddb,Vddc,Vpda,Vpdb,Vppa,Vppb,c,m1,m2,kx,ky)
 
 %定义lattice
 a1=a*[1/2,sqrt(3)/2,0];
 a2=a*[1/2,-sqrt(3)/2,0];
 %定义basis
 b1=a/sqrt(3)*[0,1,d];
 b2=a/sqrt(3)*[0,1,-d];
 %定义两个比例
 r=1/sqrt(1+d^2);
 v=1/sqrt(1+4*d^2);
 %定义k向量
 k=4*pi/(sqrt(3)*a)*[kx,ky,0];

 %不考虑其他效应，定义最原始的TB矩阵
 tb_no_soc=zeros(9,9);
 %定义对角元
 tb_no_soc(1,1)=e1;
 tb_no_soc(2,2)=e2;
 tb_no_soc(3,3)=e3;
 tb_no_soc(4,4)=e4;
 tb_no_soc(5,5)=e5;
 tb_no_soc(6,6)=e6;
 tb_no_soc(7,7)=e4;
 tb_no_soc(8,8)=e5;
 tb_no_soc(9,9)=e6;
 %定义其他元素
 %{
   矩阵前三行代表的尸V原子的dx2-dy2、dxy、dz2轨道
   对角项只考虑onsite,不考虑临近V原子的hop，因为onsite占据主导地位
   非对角项不考虑Onsite，因为onsite原子的不同磁量子数点轨道正交（？），应该为0
   下面来写tb_no_soc矩阵前三行前三列
 %}
 tb_no_soc(1,2)=sqrt(3)/16*(3*Vdda-4*Vddb+Vddc)*(exp(i*dot(k,a2))+exp(-i*dot(k,a2))-exp(i*dot(k,a1))-exp(-i*dot(k,a1)));
 tb_no_soc(2,1)=conj(tb_no_soc(1,2));
 tb_no_soc(1,3)=(Vddc-Vdda)*(-1/8*(exp(i*dot(k,a2))+exp(-i*dot(k,a2))+exp(i*dot(k,a1))+exp(-i*dot(k,a1)))+1/4*(exp(-i*dot(k,a1+a2))+exp(i*dot(k,a1+a2))));
 tb_no_soc(3,1)=conj(tb_no_soc(1,3));
 tb_no_soc(2,3)=3/8*(Vddc-Vdda)*(-exp(i*dot(k,a2))-exp(-i*dot(k,a2))+exp(i*dot(k,a1))+exp(-i*dot(k,a1)));
 tb_no_soc(3,2)=conj(tb_no_soc(2,3));
 %{
    矩阵的第4-6行分别表示的是上面一层Se原子点px py pz轨道，第7-9行是下面一层Se原子的px py pz轨道
    只考虑同种原子的相互作用
    同样，对角项只考虑onsite，非对角项不考虑onsite
    以下就是
 %}
 tb_no_soc(4,5)=-sqrt(3)/4*(Vppb-Vppa)*(-exp(i*dot(k,a2))-exp(-i*dot(k,a2))+exp(i*dot(k,a1))+exp(-i*dot(k,a1)));
 tb_no_soc(5,4)=conj(tb_no_soc(4,5));
 tb_no_soc(4,6)=0;
 tb_no_soc(6,4)=0;
 tb_no_soc(5,6)=0;
 tb_no_soc(6,5)=0;
 tb_no_soc(7,8)=tb_no_soc(4,5);
 tb_no_soc(8,7)=conj(tb_no_soc(7,8));
 tb_no_soc(7,9)=0;
 tb_no_soc(8,9)=0;
 tb_no_soc(9,7)=0;
 tb_no_soc(9,8)=0;
 %{
   下面几项都是d轨道与p轨道的hop
 %}
 tb_no_soc(1,4)=(3/8*r^3*Vpda+sqrt(3)/2*r*(1-r^2/2)*Vpdb)*(exp(i*dot(k,a2+b1))-exp(i*dot(k,b1-a1)));
 tb_no_soc(4,1)=conj(tb_no_soc(1,4));
 tb_no_soc(1,5)=-(sqrt(3)/2*r^3*Vpda+r*(1-r^2)*Vpdb)*exp(i*dot(k,b1))-(sqrt(3)/8*r^3*Vpda-r/2*(1+r^2/2)*Vpdb)*(exp(i*dot(k,a2+b1))+exp(i*dot(k,b1-a1)));
 tb_no_soc(5,1)=conj(tb_no_soc(1,5));
 tb_no_soc(1,6)=(sqrt(3)/4*d*r^3*Vpda-d/2*r^3*Vpdb)*(exp(i*dot(k,a2+b1))+exp(i*dot(k,b1-a1)))+(-sqrt(3)/2*d*r^3*Vpda+d*r^3*Vpdb)*exp(i*dot(k,b1));
 tb_no_soc(6,1)=conj(tb_no_soc(1,6));
 tb_no_soc(1,7)=tb_no_soc(1,4);
 tb_no_soc(7,1)=conj(tb_no_soc(1,7));
 tb_no_soc(1,8)=tb_no_soc(1,5);
 tb_no_soc(8,1)=conj(tb_no_soc(1,8));
 tb_no_soc(1,9)=-tb_no_soc(1,6);
 tb_no_soc(9,1)=conj(tb_no_soc(1,9));
 tb_no_soc(2,4)=r*Vpdb*exp(i*dot(k,b1))+(-3*sqrt(3)/8*r^3*Vpda-r/2*(1-3*r^2/2)*Vpdb)*(exp(i*dot(k,a2+b1))+exp(i*dot(k,b1-a1)));
 tb_no_soc(4,2)=conj(tb_no_soc(2,4));
 tb_no_soc(2,5)=(3/8*r^3*Vpda+sqrt(3)/2*r*(1-r^2/2)*Vpdb)*(exp(i*dot(k,a2+b1))-exp(i*dot(k,b1-a1)));
 tb_no_soc(5,2)=conj(tb_no_soc(2,5));
 tb_no_soc(2,6)=-sqrt(3)/4*d*r^3*(sqrt(3)*Vpda-2*Vpdb)*(exp(i*dot(k,a2+b1))-exp(i*dot(k,b1-a1)));
 tb_no_soc(6,2)=conj(tb_no_soc(2,6));
 tb_no_soc(2,7)=tb_no_soc(2,4);
 tb_no_soc(7,2)=conj(tb_no_soc(2,7));
 tb_no_soc(2,8)=tb_no_soc(2,5);
 tb_no_soc(8,2)=conj(tb_no_soc(2,8));
 tb_no_soc(2,9)=-tb_no_soc(2,6);
 tb_no_soc(9,2)=conj(tb_no_soc(2,9));
 tb_no_soc(3,4)=-sqrt(3)*r^3/2*(sqrt(3)*d^2*Vpdb+(1-2*d^2)/2*Vpda)*(exp(i*dot(k,a2+b1))-exp(i*dot(k,b1-a1)));
 tb_no_soc(4,3)=conj(tb_no_soc(3,4));
 tb_no_soc(3,5)=r^3/2*(sqrt(3)*d^2*Vpdb+(1-2*d^2)/2*Vpda)*(exp(i*dot(k,a2+b1))+exp(i*dot(k,b1-a1))-2*exp(i*dot(k,b1)));
 tb_no_soc(5,3)=conj(tb_no_soc(3,5));
 tb_no_soc(3,6)=d*r^3*(sqrt(3)*Vpdb-(1-2*d^2)/2*Vpda)*(exp(i*dot(k,a2+b1))+exp(i*dot(k,b1-a1))+exp(i*dot(k,b1)));
 tb_no_soc(6,3)=conj(tb_no_soc(3,6));
 tb_no_soc(3,7)=tb_no_soc(3,4);
 tb_no_soc(7,3)=conj(tb_no_soc(3,7));
 tb_no_soc(3,8)=tb_no_soc(3,5);
 tb_no_soc(8,3)=conj(tb_no_soc(3,8));
 tb_no_soc(3,9)=-tb_no_soc(3,6);
 tb_no_soc(9,3)=conj(tb_no_soc(3,9));
 %{
    下面几项都是p轨道与p轨道点hop
 %}
 tb_no_soc(4,7)=Vppb;
 tb_no_soc(4,8)=0;
 tb_no_soc(4,9)=0;
 tb_no_soc(5,7)=0;
 tb_no_soc(5,8)=Vppb;
 tb_no_soc(5,9)=0;
 tb_no_soc(6,7)=0;
 tb_no_soc(6,8)=0;
 tb_no_soc(6,9)=Vppa;
 tb_no_soc(7,4)=conj(tb_no_soc(4,7));
 tb_no_soc(8,4)=conj(tb_no_soc(4,8));
 tb_no_soc(9,4)=conj(tb_no_soc(4,9));
 tb_no_soc(7,5)=conj(tb_no_soc(5,7));
 tb_no_soc(8,5)=conj(tb_no_soc(5,8));
 tb_no_soc(9,5)=conj(tb_no_soc(5,9));
 tb_no_soc(7,6)=conj(tb_no_soc(6,7));
 tb_no_soc(8,6)=conj(tb_no_soc(6,8));
 tb_no_soc(9,6)=conj(tb_no_soc(6,9));

 %定义SOC矩阵
 soc=zeros(18,18);
 soc(1,2)=-2*i;
 soc(2,1)=conj(soc(1,2));
 soc(10,11)=2*i;
 soc(11,10)=conj(soc(10,11));
 
 soc(4,5)=-i;
 soc(5,4)=conj(soc(4,5));
 soc(13,14)=i;
 soc(14,13)=conj(soc(13,14));
 soc(4,15)=1;
 soc(15,4)=conj(soc(4,15));
 soc(13,6)=-1;
 soc(6,13)=conj(soc(13,6));
 soc(5,15)=-i;
 soc(15,5)=conj(soc(5,15));
 soc(14,6)=-i;
 soc(6,14)=conj(soc(14,6));
 
 soc(7,8)=-i;
 soc(8,7)=conj(soc(7,8));
 soc(16,17)=i;
 soc(17,16)=conj(soc(16,17));
 soc(7,18)=1;
 soc(18,7)=conj(soc(7,18));
 soc(16,9)=-1;
 soc(9,16)=conj(soc(16,9));
 soc(8,18)=-i;
 soc(18,8)=conj(soc(8,18));
 soc(17,9)=-i;
 soc(9,17)=conj(soc(17,9));
 
 soc=c*soc;

 %磁效应
 mag=[m1*eye(3),zeros(3,6);zeros(6,3),m2*eye(6)];

 %定义总的矩阵
 H=[tb_no_soc+mag,zeros(9,9);zeros(9,9),tb_no_soc-mag]+soc;
end 



