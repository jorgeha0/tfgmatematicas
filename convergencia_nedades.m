syms eta sigma s real
d0=1;
b0=380;
N=@(x) 1;
g=@(x) 100+0.*x;
d= @(x) d0.*(atan(20.*(x-0.5))+pi/2);
b= @(x) b0.*(sqrt(x).*exp(-2.*x)+0.01);
m=@(x) (10.*(1-x).^2).^-1;
Gamma = @(x) exp(-int(g,s,0,x));
f=@(x) exp(-x/10/(1-x))*exp(-100*(x-eta));

inner_integral = int(f, sigma, eta, 1);

outer_integral = int(inner_integral * (sqrt(eta) * exp(-2 * eta) + 1 / 100), eta, 0, 1);
nmin=60;
nmax=65;
arrayr0=zeros(1,nmax-nmin);
a1=0; a2=1;
for n= nmin:nmax
    arraybeta=promedio(b,n,a1,a2);
    arraygamma=promedio(g,n,a1,a2);
    arrayd=promedio(d,n,a1,a2); 
    arrayr0(n-nmin+1)=r0discretoindep(n,arraygamma,arrayd,arraybeta);
end
plot(linspace(nmin,nmax,nmax-nmin+1),arrayr0,'.');

function promedio=promedio(f,n,a1,a2)
    arrayx=linspace(a1,a2,n);
    arraydiscreto=zeros(1,n);
    for i=1:(length(arrayx)-1)
        integral=quad(f,arrayx(i),arrayx(i+1));
        arraydiscreto(i)=integral/((a2-a1)/n);
    end
    arraydiscreto(n)=quad(f,arrayx(n-1),a2)/((a2-a1)/n);
    promedio=arraydiscreto;
end
