AC_INIT([monoreg], 2.0)
# The following is copied from configure.ac of package `gsl': 
# Following lines kindly supplied by Dirk Eddelbuettel
## Use gsl-config to find arguments for compiler and linker flags
##
## Check for non-standard programs: gsl-config(1)
AC_PATH_PROG([GSL_CONFIG], [gsl-config])
## If gsl-config was found, let's use it
if test "${GSL_CONFIG}" != ""; then
   # Use gsl-config for header and linker arguments
   GSL_CFLAGS=`${GSL_CONFIG} --cflags`
   GSL_LIBS=`${GSL_CONFIG} --libs`
   CFLAGS="$CFLAGS $GSL_CFLAGS"  # kindly supplied by Ray Brownrigg   
else
   AC_MSG_ERROR([gsl-config not found, is GSL installed?])
fi
AC_SUBST(GSL_CFLAGS) 
AC_SUBST(GSL_LIBS)
AC_CONFIG_FILES([src/Makevars])
AC_OUTPUT