def = legacy_code('initialize');
def.SFunctionName = 'umpcS';

% Is column major
def.OutputFcnSpec = 'void umpcS(single y1[3], single y2[6], single u1[3], single u2[3][3], single u3[6], single u4[3], single u5[3], single u6, single u7, single u8, single u9, single u10, single u11, single u12, single u13, single u14, single u15, single u16, single u17, single u18[3], single u19[4], single u20[4], single u21[4], single u22[6], single u23, int32 u24, single u25[90])';

def.HeaderFiles = {'uprightmpc2.h'};
def.SourceFiles = {'uprightmpc2.c', 'eigenc.c', 'funapprox.c', ...
	'auxil.c', 'error.c', ...
	'kkt.c', 'lin_alg.c', 'osqp.c', 'proj.c', 'qdldl.c', ...
	'qdldl_interface.c', 'scaling.c', 'util.c', 'workspace.c'};
% Flat structure
def.SrcPaths = {'.'};
def.IncPaths = {'.'};
def.SampleTime = 'parameterized';

def.Options.language = 'C';

%

legacy_code('sfcn_cmex_generate', def);

%

%legacy_code('compile', def);
%legacy_code('slblock_generate', def);