".onLoad" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = "MNP"))
  title <- packageDescription("MNP", lib = mylib)$Title
  ver <- packageDescription("MNP", lib = mylib)$Version
  cat(title, "\nVersion", ver, "\n")
}

