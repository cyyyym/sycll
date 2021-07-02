function [R] = Rotation(P, W, K)
TO_RAD = pi/180;
P = P*TO_RAD;
W = W*TO_RAD;
K = K*TO_RAD;
a1 = cos(P)*cos(K)-sin(P)*sin(W)*sin(K);
a2 = -cos(P)*sin(K)-sin(P)*sin(W)*cos(K);
a3 = -sin(P)*cos(W); 
b1 = cos(W)*sin(K);
b2 = cos(W)*cos(K);
b3 = -sin(W);
c1 = sin(P)*cos(K)+cos(P)*sin(W)*sin(K);
c2 = -sin(P)*sin(K)+cos(P)*sin(W)*cos(K);
c3 = cos(P)*cos(W);
R = [a1 a2 a3;b1 b2 b3;c1 c2 c3];
clear all;
clc; 
%输入控制点坐标
x=[-86.15,-53.40,-14.78,10.46]/1000; 
y=[-68.99,82.21,-76.63,64.43]/1000; 
X=[36589.41,37631.08,39100.97,40426.54]; 
Y=[25273.32,31324.51,24934.98,30319.81]; 
Z=[2195.17,728.96,2386.50,757.31]; 
%输入焦距f，外方位元素以及内方位元素初始值，n为迭代次数
x0=0.0;
y0=0.0;
phi=0.0;
omiga=0.0;
k=0.0;
m=5000.00;
f=153.24/1000;
X0=mean(X);
Y0=mean(Y);
Z0=mean(Z)+m*f;
%定义最小二乘所需变量；
XG=zeros(6,1);
A=zeros(8,6);
L=zeros(8,1);
n=0;
phi=phi*pi/180;
omiga=omiga*pi/180;
k=k*pi/180;
n=n+1;
%计算旋转矩阵R
a1=cos(phi)*cos(k)-sin(phi)*sin(omiga)*sin(k);
a2=-cos(phi)*sin(k)-sin(phi)*sin(omiga)*cos(k);
a3=-sin(phi)*cos(omiga);
b1=cos(omiga)*sin(k);
b2=cos(omiga)*cos(k);
b3=-sin(omiga);
c1=sin(phi)*cos(k)+cos(phi)*sin(omiga)*sin(k);
c2=-sin(phi)*sin(k)+cos(phi)*sin(omiga)*cos(k);
c3=cos(phi)*cos(omiga);
R=[a1 a2 a3;b1 b2 b3;c1 c2 c3];
%求取最小二乘中的系数矩阵内各个值以及L矩阵的值
for i=1:1:4
    j=2*i-1;
    Z_Ava=a3*(X(1,i)-X0)+b3*(Y(1,i)-Y0)+c3*(Z(1,i)-Z0);
    A(j,1)=(a1*f+a3*x(1,i))/Z_Ava;
    A(j,2)=(b1*f+b3*x(1,i))/Z_Ava;
    A(j,3)=(c1*f+c3*x(1,i))/Z_Ava;
    A(j+1,1)=(a2*f+a3*y(1,i))/Z_Ava;
    A(j+1,2)=(b2*f+b3*y(1,i))/Z_Ava;
    A(j+1,3)=(c2*f+c3*y(1,i))/Z_Ava;
    A(j,4)=y(1,i)*sin(omiga)-(x(1,i)/f*(x(1,i)*cos(k)-y(1,i)*sin(k))+f*cos(k))*cos(omiga);
    A(j,5)=-f*sin(k)-x(1,i)/f*(x(1,i)*sin(k)+y(1,i)*cos(k));
    A(j,6)=y(1,i);
    A(j+1,4)=-x(1,i)*sin(omiga)-(y(1,i)/f*(x(1,i)*cos(k)-y(1,i)*sin(k))-f*sin(k))*cos(omiga);
    A(j+1,5)=-f*cos(k)-y(1,i)/f*(x(1,i)*sin(k)+y(1,i)*cos(k));
    A(j+1,6)=-x(1,i);
    
    L(j,1)=x(1,i)-(x0-f*(a1*(X(1,i)-X0)+b1*(Y(1,i)-Y0)+c1*(Z(1,i)-Z0))/Z_Ava);
    L(j+1,1)=y(1,i)-(y0-f*(a2*(X(1,i)-X0)+b2*(Y(1,i)-Y0)+c2*(Z(1,i)-Z0))/Z_Ava);
end;
%根据最小得到的公式求取观测值
XG=(inv(A'*A))*(A'*L);
  %求取地面点坐标
X0=X0+XG(1,1);
Y0=Y0+XG(2,1);
Z0=Z0+XG(3,1);
phi=phi+XG(4,1);
omiga=omiga+XG(5,1);
k=k+XG(6,1);、
%对计算误差进行判断，在误差范围内，则继续迭代，不在误差范围内，则推出循环。
while(XG(4,1)>=6.0/206265.0||XG(5,1)>=6.0/206265.0||XG(6,1)>=6.0/206265.0)
n=n+1;
a1=cos(phi)*cos(k)-sin(phi)*sin(omiga)*sin(k);
a2=-cos(phi)*sin(k)-sin(phi)*sin(omiga)*cos(k);
a3=-sin(phi)*cos(omiga);
b1=cos(omiga)*sin(k);
b2=cos(omiga)*cos(k);
b3=-sin(omiga);
c1=sin(phi)*cos(k)+cos(phi)*sin(omiga)*sin(k);
c2=-sin(phi)*sin(k)+cos(phi)*sin(omiga)*cos(k);
c3=cos(phi)*cos(omiga);
R=[a1 a2 a3;b1 b2 b3;c1 c2 c3];
for i=1:1:4
    j=2*i-1;
    Z_Ava=a3*(X(1,i)-X0)+b3*(Y(1,i)-Y0)+c3*(Z(1,i)-Z0);
    A(j,1)=(a1*f+a3*x(1,i))/Z_Ava;
    A(j,2)=(b1*f+b3*x(1,i))/Z_Ava;
    A(j,3)=(c1*f+c3*x(1,i))/Z_Ava;
    A(j+1,1)=(a2*f+a3*y(1,i))/Z_Ava;
    A(j+1,2)=(b2*f+b3*y(1,i))/Z_Ava;
    A(j+1,3)=(c2*f+c3*y(1,i))/Z_Ava;
    A(j,4)=y(1,i)*sin(omiga)-(x(1,i)/f*(x(1,i)*cos(k)-y(1,i)*sin(k))+f*cos(k))*cos(omiga);
    A(j,5)=-f*sin(k)-x(1,i)/f*(x(1,i)*sin(k)+y(1,i)*cos(k));
    A(j,6)=y(1,i);
    A(j+1,4)=-x(1,i)*sin(omiga)-(y(1,i)/f*(x(1,i)*cos(k)-y(1,i)*sin(k))-f*sin(k))*cos(omiga);
    A(j+1,5)=-f*cos(k)-y(1,i)/f*(x(1,i)*sin(k)+y(1,i)*cos(k));
    A(j+1,6)=-x(1,i);
    
    L(j,1)=x(1,i)-(x0-f*(a1*(X(1,i)-X0)+b1*(Y(1,i)-Y0)+c1*(Z(1,i)-Z0))/Z_Ava);
    L(j+1,1)=y(1,i)-(y0-f*(a2*(X(1,i)-X0)+b2*(Y(1,i)-Y0)+c2*(Z(1,i)-Z0))/Z_Ava);
end;
%完成最终迭代，输出地面点坐标及旋转矩阵
XG=(inv(A'*A))*(A'*L);
X0=X0+XG(1,1);
Y0=Y0+XG(2,1);
Z0=Z0+XG(3,1);
phi=phi+XG(4,1);
omiga=omiga+XG(5,1);
k=k+XG(6,1);
end;
format short g
R
format long g
X0,Y0,Z0
