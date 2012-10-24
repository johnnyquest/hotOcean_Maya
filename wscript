# vim: set filetype=python


def options(opt):
	opt.load("compiler_cxx")
	opt.add_option('--maya', action='store', default='2011', help='Maya version to build for')
	opt.add_option('--arch', action='store', default='x64', help='Processor architecture to build for')



def configure(cfg):
	import sys
	print sys.platform
	import os

	cwd = os.getcwd().replace('\\', '/') # FUCK windows
	is_linux = not sys.platform[0:3] == "win"

	opt = cfg.options
	arch = str(opt.arch).lower()
	
	if not is_linux:
		cfg.env['MSVC_VERSIONS'] = ['msvc 8.0', 'msvc 9.0', 'msvc 10.0']
		if arch == "x86":
			cfg.env['MSVC_TARGETS'] = ['x86']
		else:
			cfg.env['MSVC_TARGETS'] = ['x64']
		cfg.load('msvc')
	else:
		#cfg.env['CC_VERSION'] = ('4', '1', '2')
		#print cfg.env
		cfg.load('compiler_cxx')
		print cfg.env['CC_VERSION']


	env = cfg.env
	env.arch = arch
	
	if is_linux:
		env.append_value('DEFINES', 'UNIX LINUX UNIX64 LINUX_64 FUNCPROTO _GNU_SOURCE NDEBUG REQUIRE_IOSTREAM _BOOL'.split())
	else:
		if arch == "x64":
			env.append_value('DEFINES', 'WIN32 NDEBUG _WINDOWS _AFXDLL _MBCS NT_PLUGIN REQUIRE_IOSTREAM Bits64_ _WIN32 WIN32_LEAN_AND_MEAN'.split())	
		else:
			env.append_value('DEFINES', 'WIN32 NDEBUG _WINDOWS _AFXDLL _MBCS NT_PLUGIN REQUIRE_IOSTREAM _WIN32 WIN32_LEAN_AND_MEAN'.split())	
	
	includes = []

	env.maya_version = str(opt.maya).lower()
	maya_version = env.maya_version
	
	shared_path = ""
	if is_linux:
		shared_path = "/U/development/_shared/libs/"
	else:
		shared_path = "U:/development/_shared/libs/"

	includes.append(shared_path+"maya/"+maya_version+"/include/")

	# hot #1
	#
	if is_linux:
		includes.append('./3rdparty/linux/include')
	else:
		includes.append('./3rdparty/win64')

	includes.append('./3rdparty/include')


	env.append_value('INCLUDES', includes)
	
	libpath = []

	lib_arch = ""
	
	if is_linux:
		lib_arch = "lib_linux64/"
	else:
		if arch == "x64":
			lib_arch = "lib_win64/"
		else:
			lib_arch = "lib_win32/"

	libpath.append(shared_path+"maya/"+maya_version+"/"+lib_arch)

	# hot #3
	# TODO: generate these absolute paths instead
	#
	if is_linux:
		#libpath.append('./3rdparty/linux/lib/')
		libpath.append('%s/3rdparty/linux/lib/' % cwd)
	else:
		print "WHADDAFUCK:", cwd
		libpath.append('%s/3rdparty/win64/' % cwd)


	env.append_value("LIBPATH", libpath)

	libs = "Foundation OpenMaya OpenMayaUI OpenMayaAnim OpenMayaFX OpenMayaRender".split()
	libs.append("Image")

	# hot #2
	#
	if is_linux:
		libs += 'fftw3f blitz'.split()
	else:
		libs += 'libfftw3f-3 blitz'.split()
	
	#if is_linux:
	#	libs += "GL z".split()
	#else:
	#	libs += "gdi32 OpenGL32 zlib1".split()
	
	env.append_value('LIB', libs)

	
	cxxflags = []
	
	if is_linux:
		cxxflags += "-O3 -shared -m64 -fPIC -export-dynamic -Bsymbolic".split()
		cxxflags += "-fno-strict-aliasing -Wall -Wno-unused-variable -Wno-sign-compare".split()
		cxxflags += "-fopenmp".split()
	else:
		cxxflags += '/W3 /Ox /Ob2 /Oi /Ot /D "_WINDLL" /fp:fast /MD /EHsc'.split()
	
	env.append_value('CXXFLAGS', cxxflags)

	if is_linux:
		env['cxxshlib_PATTERN']='%s.so'
	else:
		env['cxxshlib_PATTERN']='%s.mll'

	
	linkflags = []
	
	if is_linux:
		pass
	else:
		linkflags.append("/export:initializePlugin")
		linkflags.append("/export:uninitializePlugin")
		linkflags.append("/subsystem:windows")
		if arch == "x64":
			linkflags.append("/MACHINE:X64")

	env.append_value('LINKFLAGS', linkflags)



def build(bld):
	sources = []
	env = bld.env	
	is_linux = env['DEST_OS'] == 'linux'

	print env
	#print env['CC_VERSION']
	
	import os
	for top, dirs, files in os.walk('./source/deformer'):
		for nm in files:       
			fp = os.path.join(top, nm)
			spl = fp.split(".")
			if spl[-1] == "cpp":
				sources.append(fp)
	bld.shlib(source = sources, target="hotOceanDeformer")
	pass



def install(ins):
	pass

