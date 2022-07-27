
#ifdef HAVE_CONFIG_H
#include "../../config.h"
#endif

#include "./dSFMT.h"

void
get_size_of_dsfmt_t_(int* size)
{
  *size = sizeof(dsfmt_t);
}

void
dsfmt_init_gen_rand_(void* dsfmt, int* seed)
{
  dsfmt_init_gen_rand((dsfmt_t*)dsfmt, (uint32_t)*seed);
}

void
dsfmt_genrand_close1_open2_(void* dsfmt, double* value)
{
  *value = dsfmt_genrand_close1_open2((dsfmt_t*)dsfmt);
}

void
dsfmt_genrand_close0_open1_(void* dsfmt, double* value)
{
  *value = (dsfmt_genrand_close1_open2((dsfmt_t*)dsfmt) - 1.0);
}
