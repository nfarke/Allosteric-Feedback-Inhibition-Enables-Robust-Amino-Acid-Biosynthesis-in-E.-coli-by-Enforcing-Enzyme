function C=exp_colormap(colors, N)
% a simple function that returns a colormap, C, for visualizing gene
% expression.  C is just a N x 3 matrix [R G B] describing the range of color values.
%
% example usage:
%    C = exp_colormap('blue-yellow',64);
%    colormap(C);
% 
% called without any arguments, it returns a [3 x 64] green-red colormap.
%
% options for 'colors' are:
%
%  'blue-yellow'
%  'green-red'
%  'yellow'
%
% 'N' represents the number of degrees of color to use.  the default is 64.
%
% generally speaking, the two-colored maps are appropriate for visualizing
% expression data which is normalized to a reference condition, e.g. to show
% whether expression is lower (blue or green) or higher (yellow or red) than
% the reference.
%
% the single-color yellow map ('yellow') is appropriate for displaying
% levels of gene expression which are not compared (necessarily) to a single
% reference, and this is similar to the colormap used in the D-chip
% expression analysis software.
% 

% the colormaps returned range monotonically.
%  
  
if 1 ~= exist('colors') 
  colors = 'green-red';
end

if 1 ~= exist('N') 
  N = 64;
end

X = [0.5: -1/(N-1):-0.5];
X = abs(X).*2;

switch colors
case {'green-red'}
  R = [X(:, 1:N/2) zeros(1,N/2)];
  G = [zeros(1, N/2) X(:,(N/2 + 1):N)];
  B = zeros(1,N);

case {'blue-yellow'} 
  R = [zeros(1,N/2) X(:,(N/2 + 1):N)];
  B = [X(:,1:N/2) zeros(1,N/2)];
  G = [zeros(1,N/2) X(:,(N/2 + 1):N)];
  
    case {'blue-orange'}
 R = [zeros(1,N/2) X(:,(N/2 + 1):N)];
  B = [X(:,1:N/2) zeros(1,N/2)];
  G = [X(:,1:N/2)/1.3 X(:,(N/2 + 1):N)/1.8];

  case {'purple-yellow'} 
  R = [X];
  B = [X(:,1:N/2) zeros(1,N/2)];
  G = [zeros(1,N/2) X(:,(N/2 + 1):N)];
  

case {'yellow'} 
  X = [0:1/(N - 1):1];
  R = X;
  B = zeros(1,N);
  G = X;

  case {'purple1'} 
  X = [0:1/(N - 1):1];
  R = X;
  B = X(:,1:N);
  G = zeros(1,N);

  
  
  
  case {'blue'} 
  
  R = [247	222	198	158	107	66	33	8	8]./255;
  G = [251	235	219	202	174	146	113	81	48]./255;
  B = [255	247	239	225	214	198	181	156	107]./255;    
  
  case {'blue2'} 
%   
%   R = ones(1,11).*0;
%   G = [0:0.085:0.85];%[0.1:0.1:1];
%   B = [0:0.1:1];    
  
  inc=1./N;

  R = ones(1,N+1).*0;
  %G = [0:(0.85.*inc):(0.85)];%[0.1:0.1:1];
  B = [0:inc:1];   
  G = B.*0.85;
   case {'purple'} 
  
  R = [252 239 218 188 158 128 106 84 63]./255;

  G = [251 237 218 189 154 125 81 39 0]./255;
  
  B = [253 245 235 220 200 186 163 143 125]./255;    
  
  
  case {'GreenBlue'} 
  
  R = [247	224	204	168	123	78	43	8	8];

  G = [252	243	235	221	204	179	140	104	64];
  
  B = [240	219	197	181	196	211	190	172	129];    
  
  
  
otherwise
 error([colors ' is not a known option for coloring']);
end

C = [R' G' B'];
