#this script unloads all libraries expect the basic ones in order to install new packages without closing the R session, sometimes still does not work.

unload_libraries <- function() {
  # Get the list of currently loaded packages
  loaded_packages <- (.packages(all.available = TRUE))
  
  # Define a vector of base packages that should not be unloaded
  base_packages <- c("base", "datasets", "graphics", "grDevices", "methods", "stats", "utils")
  
  # Unload each package one by one
  for (pkg in loaded_packages) {
    if (!(pkg %in% base_packages)) {
      pkg_detach <- paste0("package:", pkg)
      if (pkg_detach %in% search()) {
        detach(pkg_detach, unload = TRUE, character.only = TRUE)
      }
    }
  }
}

# Call the function to unload non-base libraries
unload_libraries()
