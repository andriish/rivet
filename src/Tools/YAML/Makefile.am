lib_LTLIBRARIES = libRivetYAML.la
libRivetYAML_la_SOURCES = \
  alias.cpp conversion.cpp emitterstate.cpp exp.cpp map.cpp ostream.cpp \
  parserstate.cpp scalar.cpp scanscalar.cpp sequence.cpp stream.cpp \
  content.cpp emitter.cpp emitterutils.cpp  iterator.cpp  node.cpp  \
  parser.cpp regex.cpp scanner.cpp scantoken.cpp simplekey.cpp

check_PROGRAMS = testYAML
testYAML_SOURCES = test.cc
testYAML_LDADD = -L$(builddir) -lRivetYAML
TESTS = testYAML
