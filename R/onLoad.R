".onLoad" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = "MNP"))
  title <- packageDescription("MNP", lib = mylib)$Title
  ver <- packageDescription("MNP", lib = mylib)$Version
  url <- packageDescription("MNP", lib = mylib)$URL
  cat(title, "\nVersion:", ver, "\nURL:", url, "\n")
}

