lib_LTLIBRARIES = libyaml-cpp.la
libyaml_cpp_la_SOURCES = \
  alias.cpp conversion.cpp emitterstate.cpp exp.cpp map.cpp ostream.cpp \
  parserstate.cpp scalar.cpp scanscalar.cpp sequence.cpp stream.cpp \
  content.cpp emitter.cpp emitterutils.cpp  iterator.cpp  node.cpp  \
  parser.cpp regex.cpp scanner.cpp scantoken.cpp simplekey.cpp

yamlincludedir = $(includedir)/yaml-cpp
yamlinclude_HEADERS = \
  alias.h emittermanip.h indentation.h node.h regex.h setting.h \
  content.h emitterstate.h iterator.h noncopyable.h scalar.h stlemitter.h \
  conversion.h emitterutils.h iterpriv.h ostream.h scanner.h stream.h \
  crt.h exceptions.h ltnode.h parser.h scanscalar.h token.h \
  emitter.h exp.h map.h parserstate.h sequence.h yaml.h

check_PROGRAMS = testYAML
testYAML_SOURCES = test.cc
testYAML_LDADD = -L$(top_builddir)/src/Tools/YAML -lyaml-cpp
TESTS = testYAML
