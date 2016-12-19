#if defined HAVE_CONFIG_H
#include "config.h"
#endif

module psp_zList
  use pspNode, ONLY: LIST_DATA => zNode

  include "pspLinkedlist.F90"

end module psp_zList
