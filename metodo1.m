function g = metodo1(s,met,soglia)
%compute edge-stopping function
%Input variables
%s: values to evaluate (array of length n+2)
%met: choice of the function
%     0 lorentzian function
%       g = exp(-(s.*s)./(2*sigma*sigma));
%     1 weickert function
%       g = 1 - (s.*s>0).*exp(-3.315./((s/sigma).^8));
%     2 tuckey's biweight
%       g = (abs(s)<sigma).*((1/2)*((1-(s./sigma).^2)).^2);
%soglia: edge-stopping parameter
%Output variables
%g: edge-stopping function evaluated in the given points

if met == 0 %lorenzian
    sigma=soglia;
    g = exp(-(s.*s)./(2*sigma*sigma));
end
if met==1 %weickert
    sigma=soglia;
    g = 1 - (s.*s>0).*exp(-3.315./((s/sigma).^8));
end
if met==2 %tuckey's biweight
    sigma = soglia*sqrt(5);
    g = (abs(s)<sigma).*((1/2)*((1-(s./sigma).^2)).^2);
end