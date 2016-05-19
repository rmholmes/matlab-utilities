function map = redblue_light(m);
% ctopo.m - colormap for positive/negative data with light values in between 
%
% use:  map = redblue_light(m);
% input:
%    input - colormap length (64)
%
% output:
%    output - colormap
%
% example:
%    cm1 = redblue_light(255);
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     22-May-2009
% modified: 22-May-2009
%
% obs: original from cushman-roisin book cd-rom
%     the original length was 94 instead of 255

map = [...
   127     0     0
   136     0     0
   145     0     0
   154     0     0
   163     0     0
   173     0     0
   182     0     0
   191     0     0
   200     0     0
   209     0     0
   218     0     0
   227     0     0
   236     0     0
   245     0     0
   255     0     0
   255     0     0
   255     8     8
   255    16    16
   255    24    24
   255    32    32
   255    41    41
   255    49    49
   255    57    57
   255    65    65
   255    74    74
   255    82    82
   255    90    90
   255    98    98
   255   106   106
   255   115   115
   255   123   123
   255   131   131
   255   139   139
   255   148   148
   255   156   156
   255   164   164
   255   172   172
   255   180   180
   255   189   189
   255   197   197
   255   205   205
   255   213   213
   255   222   222
   255   230   230
   255   238   238
   255   246   246
   255   255   255
   255   255   255
   246   246   255
   238   238   255
   230   230   255
   222   222   255
   213   213   255
   205   205   255
   197   197   255
   189   189   255
   180   180   255
   172   172   255
   164   164   255
   156   156   255
   148   148   255
   139   139   255
   131   131   255
   123   123   255
   115   115   255
   106   106   255
    98    98   255
    90    90   255
    82    82   255
    74    74   255
    65    65   255
    57    57   255
    49    49   255
    41    41   255
    32    32   255
    24    24   255
    16    16   255
     8     8   255
     0     0   255
     0     0   255
     0     0   245
     0     0   236
     0     0   227
     0     0   218
     0     0   209
     0     0   200
     0     0   191
     0     0   182
     0     0   173
     0     0   163
     0     0   154
     0     0   145
     0     0   136
     0     0   127];

map = map./255;

if nargin > 0
  if m ~=94
    map = interp1(1:94,map,linspace(1,94,m));
  end
end
