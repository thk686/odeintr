contents = readLines("src/RcppExports.cpp")
protos = grep("^RcppExport", contents, value = TRUE)
funcs = sub("^RcppExport SEXP (.*)\\(.*", "\\1", protos)
nargs = sapply(gregexpr("\\<SEXP\\>", protos), length) - 1
gen_call_method = function(fname, nargs)
{
  if (is.na(fname) || is.na(nargs))
    return("\t{NULL, NULL, 0}};")
  return(paste0('\t{"', fname, '", (DL_FUNC) &',
                fname, ', ', nargs, '}'))
}
funcs = c(funcs, NA); nargs = c(nargs, NA)
cms = lapply(1:length(funcs), function(i) gen_call_method(funcs[i], nargs[i]))
  
cat("#include <R.h>\n", file = "src/init.c")
cat("#include <Rinternals.h>\n", file = "src/init.c", append = TRUE)
cat("#include <R_ext/Rdynload.h>\n\n", file = "src/init.c", append = TRUE)

protos = sub(" \\{$", ";", sub("RcppExport", "extern", protos));
protos = paste0(protos, collapse = "\n")
cat(protos, file = "src/init.c", append = TRUE)

cat("\n\nR_CallMethodDef callMethods[]  = {\n", file = "src/init.c", append = TRUE)
cat(paste(cms, collapse = ",\n"), file = "src/init.c", append = TRUE)
cat("\n\nvoid\nR_init_myLib(DllInfo *info)\n{\n", file = "src/init.c", append = TRUE)
cat("\tR_registerRoutines(info, NULL, callMethods, NULL, NULL);\n", file = "src/init.c", append = TRUE)
cat("\tR_useDynamicSymbols(info, FALSE);\n", file = "src/init.c", append = TRUE)
cat("\tR_forceSymbols(info, TRUE);\n}", file = "src/init.c", append = TRUE)

