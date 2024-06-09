#ifndef SL_GRAPH_INCLUDED
#define SL_GRAPH_INCLUDED
//
// The graph structure
// Note: edges are numbered starting from 0
// IT IS IMPORTANT THAT ALL WEDGES ARE INSERTED BEFORE FIXED
//

typedef struct {
  int     nnodes;          /* number of nodes */
  int     nedges;          /* number of variables (edges) */
  unsigned short     *head;            /* head of edge e */
  unsigned short     *tail;            /* tail of edge e */
  int     *firstout;        /* first edge adjacent to node v */
  int     *nextout;         /* next edge for tail node of edge e */
  int     *firstin;         /* first edge adjacent to node v */
  int     *nextin;          /* next edge for tail node of edge e */
  char    *fixed;           /* whether arc is fixed */
  char    *unused;          /* whether arc is unused */
  int     *warc;            /* arc in wgraph (only for warcs) */
  double  *xval;            /* xval of arc (only for warcs) */
} Graph;
	
/* useful macros for manipulating the graph */

#define allocate_graph(g,n,nwedges,numfixed)                                \
  (g).head            = (unsigned short *)    malloc((nwedges+numfixed)*sizeof(unsigned short));    \
  (g).tail            = (unsigned short *)    malloc((nwedges+numfixed)*sizeof(unsigned short));    \
  (g).firstout        = (int *)    malloc((n)*sizeof(int));    \
  (g).nextout         = (int *)    malloc((nwedges+numfixed)*sizeof(int));    \
  (g).firstin         = (int *)    malloc((n)*sizeof(int));    \
  (g).nextin          = (int *)    malloc((nwedges+numfixed)*sizeof(int));    \
  (g).fixed           = (char *)   malloc((nwedges+numfixed)*sizeof(char));   \
  (g).unused          = (char *)   malloc((nwedges+numfixed)*sizeof(char));   \
  (g).warc            = (int *)    malloc((nwedges)*sizeof(int));             \
  (g).xval            = (double *) malloc((nwedges)*sizeof(double));          \

#define free_graph(g)           \
  free ( (g).head );            \
  free ( (g).tail );            \
  free ( (g).warc );            \
  free ( (g).firstout );        \
  free ( (g).nextout );         \
  free ( (g).firstin );         \
  free ( (g).nextin );          \
  free ( (g).fixed );           \
  free ( (g).xval );            \
  free ( (g).unused );          \

#define init_graph(g,n)                    \
  { int i;                                 \
    for (i=0; i<n; i++) {                  \
      (g).firstout[i] = -1;                \
      (g).firstin[i] = -1;                 \
    }                                      \
    (g).nnodes = n;                       \
    (g).nedges = 0;                       \
  }

#define insert_wedge(g,u,v,w_arc,x_val) \
  { int e = (g).nedges++;         \
    (g).head[e] = v;              \
    (g).tail[e] = u;              \
    (g).warc[e] = w_arc;          \
    (g).xval[e] = x_val;          \
    (g).fixed[e] = 0;             \
    (g).nextout[e] = (g).firstout[u];  \
    (g).firstout[u] = e;               \
    (g).nextin[e] = (g).firstin[v];    \
    (g).firstin[v] = e;                \
  }

#define insert_fixededge(g,u,v) \
  { int e = (g).nedges++;       \
    (g).head[e] = v;            \
    (g).tail[e] = u;            \
    (g).fixed[e] = 1;           \
    (g).nextout[e] = (g).firstout[u];  \
    (g).firstout[u] = e;               \
    (g).nextin[e] = (g).firstin[v];    \
    (g).firstin[v] = e;                \
  }
  

/*
#define for_neighbors(g,u,v,e)                                               \
  for (e = (g).first[u], v = ((u==(g).head[e]) ? (g).tail[e] : (g).head[e]); \
       e;                                                                    \
       e = ((u==(g).head[e]) ? (g).nexth[e] : (g).nextt[e]),                 \
       v = ((u==(g).head[e]) ? (g).tail[e] : (g).head[e]))

#define for_star(g,u,e)  \
  for (e = (g).first[u]; \
       e;                \
       e = ((u==(g).head[e]) ? (g).nexth[e] : (g).nextt[e]))

#define for_edges(g,e) \
  for (e = 0; e<*(g).nedges; e++)

#define for_pairs(g,t,h,e)                      \
  for (e = 0, t = (g).tail[e], h = (g).head[e]; \
  e < *(g).nedges;                             \
  e++, t = (g).tail[e], h = (g).head[e])

#define for_nodes(g,v) \
  for (v = 0; v < *(g).nnodes; v++)

*/
#endif