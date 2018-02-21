// calc Ibg, background Na, K, Ca 

#include "sa.h"

//////// Na ////////////
REAL
ibgna(struct State *S, struct State *Sn, REAL ht,  struct Cpar *C, struct Nernst *nernst )
{
  return C->gbna*(S->E -nernst->ENa);

} /** ibgna **/


//////// K //////////////
REAL
ibgk(struct State *S, struct State *Sn, REAL ht,  struct Cpar *C, struct Nernst *nernst )
{
  return C->gbk*(S->E -nernst->EK);

} /** ibgk **/


///////// Ca /////////////
REAL
ibgca(struct State *S, struct State *Sn, REAL ht,  struct Cpar *C, struct Caintra_state *Ca , struct Nernst *nernst )
{
  return C->gbca*(S->E -nernst->ECa);

} /** ibgca **/
