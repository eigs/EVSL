%module evsl
%{
#include "struct.h"
#include "evsl.h"
#define SWIG_FILE_WITH_INIT
%}
%include "numpy.i"
%include "struct.h"
%include "evsl.h"

%init %{
  import_array();
%}
%apply (int DIM1  , long* INPLACE_ARRAY1)
      {(int length, long* data          )};
%apply (long** ARGOUTVIEW_ARRAY1, int* DIM1  )
{(long** data , int* length)};








