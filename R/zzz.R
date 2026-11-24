# Assisted-by: Claude (Anthropic). Reviewed and validated by S. Gregoricchio.

#' @title .onUnload
#'
#' @description Releases the compiled code of the package when its namespace is unloaded. Without it a package detached with \code{unload = TRUE} and reinstalled in the same session would go on running the CRAM reader of the previous build.
#'
#' @param libpath String with the path of the library the package was loaded from.
#'
#' @author Sebastian Gregoricchio
#'
#' @noRd

.onUnload <-
    function(libpath) {
        library.dynam.unload("RegionSetDE", libpath)
    }
