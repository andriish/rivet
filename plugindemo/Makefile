# This expects all Rivet environment variables to be set up correctly already
# (use the rivetenv.{c}sh script created by the bootstrap script)

shared_flags = -shared -fPIC
# on Mac OSX you want to use these flags instead:
#shared_flags = -undefined dynamic_lookup -bundle

RivetMyAnalysis.so: MyAnalysis.cc
	g++ $(shared_flags) `rivet-config --cppflags` MyAnalysis.cc -o RivetMyAnalysis.so

all: RivetMyAnalysis.so

install: all

clean:
	rm -f RivetMyAnalysis.so
