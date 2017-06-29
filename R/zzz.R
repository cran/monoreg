
.onLoad <- function(lib, pkg) {
    library.dynam("monoreg", pkg, lib)
}

.onUnload <- function(lib) {
    library.dynam.unload("monoreg", lib)
}
