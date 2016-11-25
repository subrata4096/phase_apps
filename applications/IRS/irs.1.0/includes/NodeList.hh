#ifndef NodeList_H
#define NodeList_H

#ifdef SILO
#include <silo.h>
#endif

#include "Extents.h"

/* RGST_START */

#include "irsdefs.h"

/* RGST_END */

/*------------------------------------------------------------------------------
- DEFINES
------------------------------------------------------------------------------*/

/* flags for ndx list type generation */

#define NDX_OVERLAP          0  /* generate index list with overlaps at gblk bnds */
#define NDX_NO_OVERLAP       1  /* don't generate index list with overlaps at gblk bnds */
#define NDX_REAL             2  /* generate index list which only contains real zones/nodes */
#define NDX_USER_PHONY       3  /* include user block phonies */
#define NDX_GLBL_PHONY       4  /* include global block phonies */
#define NDX_ZONAL            5  /* zonal index list */
#define NDX_NODAL            6  /* nodal index list */

/*------------------------------------------------------------------------------
- STRUCTURE
------------------------------------------------------------------------------*/

/* RGST_START */

struct NodeWindow_t_ {

   char label[MAXLINE]; 

   int ublk;             /* user block which this list pertains to */

   int imin;             /* extents of node list */
   int imax;
   int jmin;
   int jmax;
   int kmin;
   int kmax;

   int in[4][3];    /* extents and stride info */
   int rank;             /* rank of extents */

   int len;              /* length of blk, ndx, and pos */
   int len_total;        /* length of global list */
   int logical;          /* optional logical direction */

   int *blk;             /* local block number for each ndx point */
   int *ndx;             /* local index list  */
   int *pos;             /* position within global list */
   int *udx;             /* user index list list */

   struct NodeWindow_t_ *next;

} ;

typedef struct NodeWindow_t_ NodeWindow_t;

struct NodeList_t_ {

   char name[MAXLINE];

   NodeWindow_t *list;

} ;

typedef struct NodeList_t_ NodeList_t;

/* RGST_END */


/*------------------------------------------------------------------------------
- PROTOTYPES
------------------------------------------------------------------------------*/

void NodeList_add(NodeList_t *new_ndx);
void NodeList_free(void);
void NodeList_make( NodeWindow_t *ndxin, char *name );
NodeList_t *NodeList_find( char *name );
void NodeList_addnextseq( char *name, NodeWindow_t *ndxin , char *result);
int NodeList_read( void );
void NodeList_print(NodeList_t *NodeList);
int printnodelist(void);

#ifdef SILO
int NodeList_wtsilo(DBfile *idbid);
int NodeList_rdsilo(DBfile *idbid, int ublk0);
int NodeWindow_rdsilo(DBfile *idbid, char *name, NodeWindow_t **ndxin, int ublk0);
int NodeWindow_rdsilo(DBfile *idbid, char *name, NodeWindow_t **ndxin, int ublk0);
int NodeWindow_wtsilo(DBfile *idbid, char *name, NodeWindow_t *ndxin);
#endif

void NodeList_pack( int **data, int *outlen, int *stroutlen, int oldlen );
void NodeList_cpack( char **data );
void NodeList_unpack( int *data, int len, char *ndx_names );
char *NodeList_cunpack( char *data, int len );
void NodeWindow_setext(Extents_t *ext, int ublk, int ndx_incl, int ndx_cent) ;
int  NodeWindow_read(char *name, NodeWindow_t **in_ndx) ;
void NodeWindow_copy(NodeWindow_t *ndxout, NodeWindow_t *ndxin) ;
void NodeWindow_add(NodeWindow_t *new_ndx, NodeWindow_t **list) ;
void NodeWindow_free(NodeWindow_t *ndxin) ;
void NodeWindow_freendx(NodeWindow_t *ndxin);
int  NodeWindow_getndx(NodeWindow_t *inndx, int gblk_in, int ndx_ovlp, int ndx_incl,
                   int ndx_cent );
int  NodeWindow_getlen(NodeWindow_t *inndx, int gblk_in, int ndx_ovlp, int ndx_incl,
                   int ndx_cent);
int  NodeWindow_range(NodeWindow_t *ndxin1, NodeWindow_t *ndxin2, int plusmin1, int plusmin2);
void NodeWindow_setndx(int ublk_in, int gblk_in, int i, int j, int k, int* n,
                 NodeWindow_t *inndx, int ndx_ovlp, int ndx_incl, int ndx_cent);
int NodeWindow_compact(NodeWindow_t *inndx, NodeWindow_t *out_ndx,
                    int gblk_in,
                    int ndx_ovlp,
                    int ndx_incl,
                    int ndx_cent ) ;
void NodeWindow_copyall(NodeWindow_t *in_ndx, NodeWindow_t **out_ndx);

void NodeWindow_make( NodeWindow_t **ndxin, int imin, int imax,
                   int jmin, int jmax, int kmin, int kmax, int ublk,
                   char *name)  ;


void NodeWindow_buildndx( NodeWindow_t *ndxin, int ndx_incl, int ndx_cent,
                      int gblk, int *inlen, int **inndx);



void NodeWindow_intersect(char *name, NodeWindow_t *ndxin1, 
                          NodeWindow_t *ndxin2,
                          NodeWindow_t **ndxout) ;

void NodeList_rename( char *name1, char *name2 );

void NodeList_del( char *name);

void NodeWindow_getstr(int ndx, char *string, int gblk); 

int NodeList_plnl( void ) ;
int NodeList_plnloff( void ) ;
void NodeList_plot( int iflab );
int NodeWindow_fastndx(NodeWindow_t *in_ndx, int gblk, int ndx_incl, int ndx_cent);

#endif

/*------------------------------------------------------------------------------
- GLOBAL VARS
------------------------------------------------------------------------------*/

/*******************************************************************************
* END OF FILE
*******************************************************************************/
