function [y] = epsilon (x1, x2, x, epsilon0, epsilon1)

% La razón por la que esta M-función ha sido concebida no es otra que la de
% otorgar un valor de la permitivad eléctrica (o cómo se llame, todos
% sabemos a qué nos referimos). En concreto, si estamos en la lámina, debe
% darnos el valor de la epsilon de la lámina y si estamos en el vacío, pues
% la del vacío.

d=length(x);

y=zeros(1,d);

for i=1:d

if (x(i)>x1 && x2>x(i))
    
    y(i)=epsilon1;
    
else
    
    y(i)=epsilon0;
    
end

end

end