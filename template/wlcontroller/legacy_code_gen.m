def = legacy_code('initialize');
def.SFunctionName = 'ex_wlmex';
def.OutputFcnSpec = 'void wlcWrapper(double y1[4], double u1[4], double u2[6], double u3[6], double u4[6], double u5[6], double u6[6])';

def.HeaderFiles = {'wlmex.h'};
def.SourceFiles = {'wlcontroller.cpp', 'wlqp.cpp', 'auxil.c', 'error.c', ...
	'kkt.c', 'lin_alg.c', 'osqp.c', 'proj.c', 'qdldl.c', ...
	'qdldl_interface.c', 'scaling.c', 'util.c', 'workspace.c', ...
    'wlmex.cpp'};
def.SrcPaths = {'.', 'wlqp/src/osqp'};
def.IncPaths = {'.', 'eigen', 'wlqp/include'};
def.SampleTime = 'parameterized';

def.Options.language = 'C';

%

legacy_code('sfcn_cmex_generate', def);

%

%legacy_code('compile', def);
%legacy_code('slblock_generate', def);
