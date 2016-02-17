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
        "Main 'Unix-like' part of the AO path, optionally including a /REF prefix"
        return self._basepath if keepref else self._basepath.lstrip("/REF").rstrip("/")

    def varpath(keepref=False, defaultvarid=None):
        "The basepath, plus any bracketed variation identifier"
        p = self.basepath(keepref)
        if self.varid(defaultvarid) is not None:
            p += "[%s]" % str(self.varid(defaultvarid))
        return p

    def binpath(keepref=False, defaultbinid=None, defaultvarid=None):
        "The varpath, plus any #-prefixed bin number identifier"
        p = self.varpath(keepref, defaultvarid)
        if self.binid(defaultbinid) is not None:
            p += "#%d" % self.binid(defaultbinid)
        return p

    def basepathparts(keepref=False):
        "List of basepath components, split by forward slashes"
        return self.basepath(keepref).strip("/").split("/")

    # TODO: basepathhead, basepathtail

    def dirname(keepref=False):
        "The non-final (i.e. dir-like) part of the basepath"
        return os.path.dirname(self.basepath(keepref))

    def dirnameparts(keepref=False):
        "List of dirname components, split by forward slashes"
        return self.dirname(keepref).strip("/").split("/")

    def basename():
        "The final (i.e. file-like) part of the basepath"
        return self._basename

    def varid(default=None):
        "The variation identifier (without brackets) if there is one, otherwise None"
        return self._varid if self._varid is not None else default

    def binid(default=None):
        "The bin identifier (without #) if there is one, otherwise None"
        return self._binid if self._binid is not None else default

    def isref():
        "Is there a /REF prefix in the original path?"
        return self._isref

    def istmp():
        "Do any basepath components start with an underscore, used to hide them from plotting?"
        return isTmpPath(self.basepath())
