load ijT.dat
data = ijT;

z = data(:,3);
caxis([0 1])

  for i = 1:34
    dd(1:34,i) = z( ((i-1)*34+1): (i*34) );
  endfor

pcolor(dd)

dd(1:34,35)=0;
dd(35,1:35)=0;

xlabel('x')
ylabel('y')

size(dd)

colorbar