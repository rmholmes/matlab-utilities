function map = ctopo(m);
% ctopo.m - colormap for topography
%
% use:  map = ctopo(m);
% input:
%    input - colormap length (64)
%
% output:
%    output - colormap
%
% example:
%    cm1 = ctopo(255);
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     22-May-2009
% modified: 22-May-2009
%
% obs: original matfile was only 66 length instead of 255
%

map = [...
   255   253   145
   251   251   143
   247   249   140
   243   246   138
   240   244   136
   236   242   134
   232   240   131
   228   238   129
   225   236   127
   221   234   125
   217   232   122
   213   230   120
   210   227   118
   206   225   116
   202   223   114
   198   221   111
   195   219   109
   191   217   107
   187   215   105
   183   213   102
   180   211   100
   176   208    98
   172   206    96
   168   204    93
   165   202    91
   161   200    89
   157   198    87
   153   196    84
   150   194    82
   146   192    80
   142   190    78
   138   187    76
   135   185    73
   131   183    71
   127   181    69
   124   179    67
   120   177    64
   116   175    62
   112   173    60
   109   171    58
   105   168    55
   101   166    53
    97   164    51
    94   162    49
    90   160    46
    86   158    44
    82   156    42
    79   154    40
    75   152    38
    71   149    35
    67   147    33
    64   145    31
    60   143    29
    56   141    26
    52   139    24
    49   137    22
    45   135    20
    41   133    17
    37   131    15
    34   128    13
    30   126    11
    26   124     8
    22   122     6
    19   120     4
    15   118     2
    11   116     0];

map = map./255;

if nargin > 0
  if m ~=66
    map = interp1(1:66,map,linspace(1,66,m));
  end
end