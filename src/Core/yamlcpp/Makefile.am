EXTRA_DIST = yaml-cpp-0.3.0.tar.gz

all: yaml-cpp/librivet-yaml-cpp.a
	@true

## Horrible complicated rule to force a -fPIC build, rename the library, and work around OS X sed issues
yaml-cpp/librivet-yaml-cpp.a: yaml-cpp-0.3.0.tar.gz
	test -d yaml-cpp || tar xf $^; true
	grep COMPILE_FLAGS yaml-cpp/CMakeLists.txt || $(SED) -e $$'s/LABEL_SUFFIX}"/LABEL_SUFFIX}"\\\n\\\tCOMPILE_FLAGS "-fPIC"/' \
      yaml-cpp/CMakeLists.txt > tmp; test -e tmp && mv tmp yaml-cpp/CMakeLists.txt; true
	for i in `find yaml-cpp -type f -name "*.h" -o -name "*.cpp"`; do (grep RIVET_YAML $$i > /dev/null) || \
      { $(SED) -e 's/namespace YAML/namespace RIVET_YAML/' -e 's/YAML::/RIVET_YAML::/g' $$i > tmp; \
        test -e tmp && mv tmp $$i; } done; true
	cd yaml-cpp && CC="$(CC)" CXX="$(CXX)" CXXFLAGS="$(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) $(AM_LDFLAGS) $(LDFLAGS) $(DYLDFLAGS)" $(CMAKE) . && \
      $(MAKE) && mv libyaml-cpp.a librivet-yaml-cpp.a

clean-local:
	@rm -rf yaml-cpp
