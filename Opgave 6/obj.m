function [f,df,d2f]=obj(x)

[g,dg,d2g]=obj1(x);
[h,dh,d2h]=obj2(x);

f = g-h;
df = dg-dh;
d2f = d2g-d2h;
