lib_LTLIBRARIES = libyaml-cpp.la
libyaml_cpp_la_SOURCES = \
  aliasmanager.cpp binary.cpp conversion.cpp directives.cpp emitfromevents.cpp \
  emitter.cpp emitterstate.cpp emitterutils.cpp exp.cpp iterator.cpp \
  nodebuilder.cpp node.cpp nodeownership.cpp null.cpp ostream.cpp parser.cpp \
  regex.cpp scanner.cpp scanscalar.cpp scantag.cpp scantoken.cpp simplekey.cpp \
  singledocparser.cpp stream.cpp tag.cpp

yamlincludedir = $(includedir)/yaml-cpp
yamlinclude_HEADERS = \
 aliasmanager.h anchor.h binary.h conversion.h dll.h emitfromevents.h emitter.h \
 emittermanip.h eventhandler.h exceptions.h iterator.h ltnode.h mark.h node.h \
 nodeimpl.h nodereadimpl.h nodeutil.h noncopyable.h null.h ostream.h parser.h \
 stlemitter.h stlnode.h traits.h yaml.h

dist_noinst_HEADERS = \
  collectionstack.h emitterstate.h exp.h iterpriv.h nodeownership.h ptr_vector.h \
  regeximpl.h scanscalar.h setting.h streamcharsource.h stringsource.h token.h \
  directives.h emitterutils.h indentation.h nodebuilder.h ptr_stack.h regex.h \
  scanner.h scantag.h singledocparser.h stream.h tag.h

libyaml_cpp_la_CPPFLAGS = $(AM_CPPFLAGS) -I$(builddir)/..

check_PROGRAMS = testYAML
testYAML_SOURCES = test.cc
testYAML_CPPFLAGS = -I$(builddir)/..
testYAML_LDADD = -L$(builddir) -lyaml-cpp
TESTS = testYAML
