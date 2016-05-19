function [tp,ep,A,C]=scaloa(xc,yc,x,y,t,corrlen,err)

%function [tp,ep]=scaloa(xc,yc,x,y,t,corrlen,err)

% funcao de analise objetiva para quantidade escalar
% assume funcao de correlacao espacial isotropica e gaussiana do tipo:

% *** C=(1-err)*exp(-d2/corrlen^2)***

% sendo:
% corrlen: comprimento de correlacao;
% err:     variancia do erro (aleatorio) amostral assumido % (aparece sempre como epsilon^2 nos papers)
%          (aparece sempre como epsilon^2 nos papers)
% d:       distancia radial entre os pontos de observacao

%          a rotina apenas faz a AO. Entra-se com valores de
%          corrlen e err

% tp:      vetor da quantidade gradeada
%          os dados de entrada sao t(x,y) para
%          interpolar nos pontos (xc,yc) da sua grade

% a rotina apenas faz a AO, e vc tem que entrar com valores de
% corrlen e err

% os dados de entrada sao t(x,y) para interpolar nos pontos (xc,yc) da sua grade

n=length(x);
x=reshape(x,1,n);
y=reshape(y,1,n);
t=t';

% constroe a matriz de quadrados de distancia radial entre os
% pontos de observacao

d2=((x(ones(n,1),:)'-x(ones(n,1),:)).^2+...
(y(ones(n,1),:)'-y(ones(n,1),:)).^2);
nv=length(xc);
xc=reshape(xc,1,nv);
yc=reshape(yc,1,nv);

% constroe a matriz de quadrados de distancia radial entre os
% pontos de observacao e os pontos de sua grade

dc2=((xc(ones(n,1),:)'-x(ones(nv,1),:)).^2+...
(yc(ones(n,1),:)'-y(ones(nv,1),:)).^2);

% constroe a matriz de correlacao entre as observacoes (matriz C) usando a
matriz "d"
% e a matriz de correlacao cruzada entre os pontos de observacao e
% os pontos de grade (matriz A) usando a matriz "dc"

C=(1-err)*exp(-dc2/corrlen^2);
A=(1-err)*exp(-d2/corrlen^2);

% soma a matriz diagonal associada ao erro amostral assumido
% (eh diagonal pois o erro eh assumido aleatorio, e portanto, soh se
correlaciona
% com ele na mesma localidade)

A=A+err*eye(size(A));

% Bingo...aplica-se o teorema de Gauss-Markov criando uma "funcao-peso"
% que minimiza variancia..interpolacao otima!

% Esta funcao peso eh "C*inv(A) ou em matlabe^s, escrita da forma abaixo.

if(~isempty(t))
tp=(C*(A\t))';
end

% ep eh o vetor contendo o erro medio quadratico normalizado. Se extrair
% a raiz quadrada, voce tem o erro de interpolacao percentualmente
% (obviamente).

% basta dar um "reshape" em tp e ep e vc tem os mapas da quantidade "t"
% otimamente mapeada - a matriz tp-, e do REMQ normalizado.

% note que o mapa de erro independe do valor da quantidade interpolada
% depende apenas das distancias entre os pontos de grade e os observados,
% alem, claro, do comprimento de correlacao e variancia de erro amostral
% escolhidos

if(nargout==4)
ep=1-sum(C'.*(A\C'))/(1-err);
end

% a scaloa assume que vc jah determinou estes parametros ou simplesmente os
% escolheu arbitrariamente

% caso nao tenha funcao de correlacao amostral, os tradicionais chutes iniciais sao o primeiro raio
% baroclinico local para o comp. de correlacao e um valor entre 0.1 e 0.2
% para o std do erro amostral.