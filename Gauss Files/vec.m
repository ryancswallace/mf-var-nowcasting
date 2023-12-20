function A=vec(B)

n=rows(B);
m=cols(B);
A=reshape(B,n*m,1);
