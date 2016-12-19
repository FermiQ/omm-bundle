#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module psp_dList
  use pspNode, ONLY: LIST_DATA => dNode

  include "pspLinkedlist.F90"

end module psp_dList
