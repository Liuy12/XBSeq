# ==========================================================================
# package initialization
# ==========================================================================
.onAttach = function(libname, pkgname) {
    msg = "Welcome to 'XBSeq'."
    msg = strwrap(msg, exdent=4, indent=4)
    packageStartupMessage(paste(msg, collapse="\n"), appendLF=TRUE)
}

