#if !defined(__TAO_VERSION_H)
#define __TAO_VERSION_H

/* ========================================================================== */
/* 
   Current TAO version number and release date, also listed in
    docs/changes.html
    docs/tex/manual/manual.tex 
    docs/tex/manual/intro.tex 
    docs/tex/manual/manual_tex.tex
*/
#define TAO_VERSION_NUMBER "TAO Version 1.10.1, Released May 24, 2010"
#define TAO_VERSION_RELEASE  1
#define TAO_VERSION_MAJOR    1
#define TAO_VERSION_MINOR    10
#define TAO_VERSION_SUBMINOR 1
#define TAO_PATCH_LEVEL      3
#define TAO_VERSION_(MAJOR,MINOR,SUBMINOR) \
    ((TAO_VERSION_MAJOR == (MAJOR)) &&       \
    (TAO_VERSION_MINOR == (MINOR)) &&       \
     (TAO_VERSION_SUBMINOR == (SUBMINOR)))
#define TAO_VERSION_DATE     "Jul 7, 2010"
#define TAO_AUTHOR_INFO      "The TAO Team:\
 Lois Curfman McInnes, Todd Munson, Jorge More', Jason Sarich\n\
Bug reports, questions: tao-comments@mcs.anl.gov\n\
Web page: http://www.mcs.anl.gov/tao/\n"

#endif
