def = legacy_code('initialize');
def.SFunctionName = 'ex_wlmex';

def.OutputFcnSpec = 'void wlControllerUpdate(single y1[4], single u1[4], single u2[6], single u3[6], single u4[90], single u5)';

def.HeaderFiles = {'wlcontroller.h'};
def.SourceFiles = {'wlcontroller.c', 'wlqp.c', 'eigenc.c', 'funapprox.c', ...
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
