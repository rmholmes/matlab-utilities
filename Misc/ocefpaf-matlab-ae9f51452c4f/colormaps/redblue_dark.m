function map = redblue_dark(m);
% redblue_dark - colormap for positive/negative data with dark values in between
%
% use:  map = redblue_dark(m);
% input:
%    input - colormap length (64)
%
% output:
%    output - colormap
%
% example:
%    cm1 = redblue_dark(255);
%
% author:   Filipe P. A. Fernandes
% e-mail:   ocefpaf@gmail.com
% web:      http://ocefpaf.tiddlyspot.com/
% date:     22-May-2009
% modified: 22-May-2009
%
% obs: original length 254 instead of 255
%

map = [...
     0     0   255
     0     0   252
     0     0   250
     0     0   248
     0     0   246
     0     0   244
     0     0   242
     0     0   240
     0     0   238
     0     0   236
     0     0   234
     0     0   232
     0     0   230
     0     0   228
     0     0   226
     0     0   224
     0     0   222
     0     0   220
     0     0   218
     0     0   216
     0     0   214
     0     0   212
     0     0   210
     0     0   208
     0     0   206
     0     0   204
     0     0   202
     0     0   200
     0     0   198
     0     0   196
     0     0   194
     0     0   192
     0     0   190
     0     0   188
     0     0   186
     0     0   184
     0     0   182
     0     0   180
     0     0   178
     0     0   176
     0     0   174
     0     0   172
     0     0   170
     0     0   167
     0     0   165
     0     0   163
     0     0   161
     0     0   159
     0     0   157
     0     0   155
     0     0   153
     0     0   151
     0     0   149
     0     0   147
     0     0   145
     0     0   143
     0     0   141
     0     0   139
     0     0   137
     0     0   135
     0     0   133
     0     0   131
     0     0   129
     0     0   127
     0     0   125
     0     0   123
     0     0   121
     0     0   119
     0     0   117
     0     0   115
     0     0   113
     0     0   111
     0     0   109
     0     0   107
     0     0   105
     0     0   103
     0     0   101
     0     0    99
     0     0    97
     0     0    95
     0     0    93
     0     0    91
     0     0    89
     0     0    87
     0     0    85
     0     0    82
     0     0    80
     0     0    78
     0     0    76
     0     0    74
     0     0    72
     0     0    70
     0     0    68
     0     0    66
     0     0    64
     0     0    62
     0     0    60
     0     0    58
     0     0    56
     0     0    54
     0     0    52
     0     0    50
     0     0    48
     0     0    46
     0     0    44
     0     0    42
     0     0    40
     0     0    38
     0     0    36
     0     0    34
     0     0    32
     0     0    30
     0     0    28
     0     0    26
     0     0    24
     0     0    22
     0     0    20
     0     0    18
     0     0    16
     0     0    14
     0     0    12
     0     0    10
     0     0     8
     0     0     6
     0     0     4
     0     0     2
     0     0     0
     0     0     0
     2     0     0
     4     0     0
     6     0     0
     8     0     0
    10     0     0
    12     0     0
    14     0     0
    16     0     0
    18     0     0
    20     0     0
    22     0     0
    24     0     0
    26     0     0
    28     0     0
    30     0     0
    32     0     0
    34     0     0
    36     0     0
    38     0     0
    40     0     0
    42     0     0
    44     0     0
    46     0     0
    48     0     0
    50     0     0
    52     0     0
    54     0     0
    56     0     0
    58     0     0
    60     0     0
    62     0     0
    64     0     0
    66     0     0
    68     0     0
    70     0     0
    72     0     0
    74     0     0
    76     0     0
    78     0     0
    80     0     0
    82     0     0
    85     0     0
    87     0     0
    89     0     0
    91     0     0
    93     0     0
    95     0     0
    97     0     0
    99     0     0
   101     0     0
   103     0     0
   105     0     0
   107     0     0
   109     0     0
   111     0     0
   113     0     0
   115     0     0
   117     0     0
   119     0     0
   121     0     0
   123     0     0
   125     0     0
   127     0     0
   129     0     0
   131     0     0
   133     0     0
   135     0     0
   137     0     0
   139     0     0
   141     0     0
   143     0     0
   145     0     0
   147     0     0
   149     0     0
   151     0     0
   153     0     0
   155     0     0
   157     0     0
   159     0     0
   161     0     0
   163     0     0
   165     0     0
   167     0     0
   170     0     0
   172     0     0
   174     0     0
   176     0     0
   178     0     0
   180     0     0
   182     0     0
   184     0     0
   186     0     0
   188     0     0
   190     0     0
   192     0     0
   194     0     0
   196     0     0
   198     0     0
   200     0     0
   202     0     0
   204     0     0
   206     0     0
   208     0     0
   210     0     0
   212     0     0
   214     0     0
   216     0     0
   218     0     0
   220     0     0
   222     0     0
   224     0     0
   226     0     0
   228     0     0
   230     0     0
   232     0     0
   234     0     0
   236     0     0
   238     0     0
   240     0     0
   242     0     0
   244     0     0
   246     0     0
   248     0     0
   250     0     0
   252     0     0
   255     0     0];

map = map./255;

if nargin > 0
  if m ~=254
    map = interp1(1:254,map,linspace(1,254,m));
  end
end
