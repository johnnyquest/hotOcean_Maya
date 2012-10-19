These were written and compiled in MS Visual Studio 2008.
Besides the standard Maya/Mental Ray libaries/includes you'll need the 
Houdini Ocean Toolkit source which can be found here: 
http://odforce.net/wiki/index.php/HoudiniOceanToolkit

My includes looked like this:

/I "C:\Program Files\Autodesk\Maya2009\devkit\mentalray\include" 
/I "E:\Projects\Programming\Ocean\HOT_src\3rdparty\win64" 
/I "E:\Projects\Programming\Ocean\HOT_src\3rdparty\include\loki\yasli" 
/I "E:\Projects\Programming\Ocean\HOT_src\3rdparty\include\loki\flex" 
/I "E:\Projects\Programming\Ocean\HOT_src\3rdparty\include\loki" 
/I "E:\Projects\Programming\Ocean\HOT_src\3rdparty\include" 

And the linker needed the following dependencies:

shader.lib 
mayabase.lib 
blitz.lib 
libfftw3f-3.lib