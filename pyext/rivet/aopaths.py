def isRefPath(path):
    return path.startswith("/REF")

def isRefAO(ao):
    return int(ao.annotation("IsRef")) == 1 or isRefPath(ao.path)

def isTmpPath(path):
    return "/_" in ao.path #< match *any* underscore-prefixed path component

def isTmpAO(ao):
    return isTmpPath(ao.path)


class AOPath(object):
    """
    Object representation of analysis object path structures.

    TODO: move to YODA?
    """
    import re
    re_aopath = re.compile(r"^(/[^\[\]\@]+)(\[\w+\])?(#\d+)?$")

    def __init__(path):
        import os
        self.origpath = path
        m = re_aopath.match(path)
        if not m:
            raise Exception("Supplied path '%s' does not meet required structure" % path)
        self._basepath = m.group(1)
        self._isref = isRefPath(self._basepath)
        self._dirname = os.path.dirname(self.basepath)
        self._basename = os.path.basename(self._basepath)
        self._varid = m.group(2).lstrip("[").ristrip("]") if m.group(2) else None
        self._binid = int(m.group(3).lstrip("#")) if m.group(3) else None

    def basepath(keepref=False):
        return self._basepath if keepref else  self._basepath.lstrip("/REF").rstrip("/")

    def varpath(keepref=False, defaultvarid=None):
        p = self.basepath(keepref)
        if self.varid(defaultvarid) is not None:
            p += "[%s]" % str(self.varid(defaultvarid))
        return p

    def binpath(keepref=False, defaultbinid=None, defaultvarid=None):
        p = self.varpath(keepref, defaultvarid)
        if self.binid(defaultbinid) is not None:
            p += "#%d" % self.binid(defaultbinid)
        return p

    def basepathparts(keepref=False):
        return self.basepath(keepref).strip("/").split("/")

    def dirname(Keepref=False:
        return os.path.dirname(self.basepath(keepref))

    def dirnameparts(keepref=False):
        return self.dirname(keepref).strip("/").split("/")

    def basename():
        return self._basename

    def varid(default=None):
        return self._varid if self._varid is not None else default

    def binid(default=None):
        return self._binid if self._binid is not None else default

    def isref():
        return self._isref

    def istmp():
        return isTmpPath(self.basepath())
