ACLOCAL_AMFLAGS = -I m4
SUBDIRS = src pyext include data bin doc

doc:
	cd doc && $(MAKE) doc

.PHONY : doc

clean-local:
	@rm -rf a.out
