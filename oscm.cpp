//
// oscm: one sided crossing minimization
//
// reads instance from stdin, writes solution to stdout
//

// This program is our contribution to the PACE 2024 challenge.

// License: GNU General Public License v3.0

// Bonn, Cologne, Heidelberg, June 9, 2024
//
// Michael Juenger (University of Cologne)
// Paul Juenger (University of Bonn)
// Petra Mutzel (Unibersity of Bonn)
// Gerhard Reinelt (University of Heidelberg)

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <sys/resource.h>
#include <assert.h>

#include "CbcInterface.hpp"
#include "CbcCutGenerator.hpp"
#include "OsiClpSolverInterface.hpp"

#include "graph.h"

#define SUBMISSION

//
// GLOBALS
//

// constants
#define false 0
#define true 1

// oscm globals
int pril;                    // output level
int nnod;                    // number of items
int wnvar;                   // number of variables in weighted graph
int num3cyclecuts = 0;       // total number of separated 3-cycle inequalities
int numlargecyclecuts = 0;   // total number of separated >3-cycle inequalities
int numfraccyclecuts = 0;    // total number of separated fractional cycle inequalities
int nummoebiuscuts = 0;      // total number of separated moebius inequalities
int max_gen_cycles = 3000;   // no more cycles than this

int numtopnodes;             // number of top nodes in input
int numbottomnodes;          // number of bottom nodes in input
int numedges;                // number of edges in input
int numcomponents;           // number of components
int nooutput = false;        // true if no output of the solution to stdout is wanted
int numlps = 0;              // number of lps solved
int numequal;                // number of node pairs v,w with equal crossing values regardless of ordering
int numfixed;                // number of fixed acyclic subgraph arc variables

double objconst;             // objective function constant
double eps = 0.0001;      // epsilon
double omeps = 1.0 - eps; // 1 - epsilon
double opeps = 1.0 + eps; // 1 + epsilon
double weightsum;            // sum of all arc weights
int numbbnodes = 0.0;     // number of branch & bound nodes

double currentCompLocalLowerbound = 0.0;

double processtimelimit = 1e6;     // time limit in seconds
double programstarttime = 0.0;     // starting time of the program execution
double starttime = 0.0;            // starting time of a portion of the program
double remainingtime = 0.0;        // remaining time for the computation
double elapsedtime = 0.0;          // time since beginning of computation
double inputtranstime = 0.0;       // time for reading input and transformation
double cbcinittime = 0.0;          // time for cbc to initialize model
double compstarttime = 0.0;        // time when a component bac computation starts
double exploitationtime = 0.0; // time spent in exploiting LP solutions
double separationtimecallback = 0.0;   // time spent in separation of callback inequalities
double separationtimemoebius = 0.0;    // time spent in separation of moebius inequalities
double separationtimeother = 0.0;      // time spent in separation of cycle inequalities
double componentcrosstabletime = 0.0;// time for constructing component crosstable
double componentgraphtime = 0.0;    // time for constructing component graph
double bactime = 0.0;              // time spent in the branch&cut phases
double misctime = 0.0;             // miscellaneous time not accounted for in detail
double totaltime = 0.0;            // total program execution time

// the weighted digraph used in lp/ip solving
unsigned short *wtail;                   // tail of an arc in an acyclic subdigraph instance
unsigned short *whead;                   // head of an arc in an acyclic subdigraph instance
int *wfirst;                  // first arc (v,w) in adjacency list of node v in an acyclic subdigraph instance
int *wnext;                   // next arc (v,w) in adjacency list of node v in an acyclic subdigraph instance
int *wfirstin;                // first in arc (v,w) in adjacency list of node v in an acyclic subdigraph instance
int *wnextin;                 // next in arc (v,w) in adjacency list of node v in an acyclic subdigraph instance

// the digraph of the fixed arcs
unsigned short *ftail;                  // tail of a fixed arc in an acyclic subdigraph instance
unsigned short *fhead;                  // head of a fixed arc in an acyclic subdigraph instance
int *ffirst;                 // first fixed arc (v,w) in adjacency list of node v
int *fnext;                  // next fixed arc (v,w) in adjacency list of node v

Graph* sepgraph;             // graph used for separation
Graph* sepgraphint;          // graph used for separation (only integer edges)

int *cind;                   // variable indices of a constraint to be added
int *graph_cind;             // variable indices (in sep. graph) of a constraint to be added
double *cval;                // coefficients of a constraint to be added

int *activity;               // best known linear ordering of items
int *rank;                   // ranks of items in best known linear ordering
double activityrank_num_crossings = 1e9; // number of crossings of solution stored in activity/rank

int *node;                   // array of all nodes in a component
int *poscct;                 // position of a node in compcrosstable
int *crosstable;             // table of crossings between edges incident to bottom node pairs
int *compcrosstable;         // table of crossings between edges incident to bottom node pairs in a component

double *obj;                 // objective function coefficients for weighted digraph used in lp/ip solving

int timeout = false;             // whether timeout has occurred

int iscutwidth = false; // whether to use moebius and initial heuristic
int dofixing = true;
int doarbitrary = true;
int separatemoebius; // whether to separate moebius inequalities, may change during bac

double lastRootLowerbound  = 10e9;
int numRootSmallObjChange = 0;

//
// OBTAIN PROCESS TIME
//
double processtime()
{
#ifdef SUBMISSION
    return 0;
#else
    struct rusage r;
    getrusage(RUSAGE_SELF, &r);
    const double utime = r.ru_utime.tv_sec + ((double) r.ru_utime.tv_usec / 1000000.0);
    return utime;
#endif
}

//
// PRINT PROCESS TIME
//
void printtime(double t)
{
    double     dhours,dminutes,dseconds;
    dseconds = t;
    dminutes = floor(dseconds/60.0);
    dhours   = floor(dminutes/60.0);
    dseconds = fmod(dseconds,60.0);
    dminutes = fmod(dminutes,60.0);
    printf("  %5.0lf:%02.0lf:%05.2lf",dhours,dminutes,dseconds);
}

//
// RANDOM GENERATOR (taken from Stanford GraphBase by Donald E. Knuth)
// edited in order to make compiler and reader happy
//

#define gb_next_rand() (*gb_fptr>=0?*gb_fptr--:gb_flip_cycle())
#define mod_diff(x,y)(((x)-(y))&0x7fffffff)
#define two_to_the_31 ((unsigned long)0x80000000)

long A[56]= {-1};
long*gb_fptr= A;

long gb_flip_cycle()
{
    long*ii,*jj;
    for(ii= &A[1],jj= &A[32];jj<=&A[55];ii++,jj++)
        *ii= mod_diff(*ii,*jj);
    for(jj= &A[1];ii<=&A[55];ii++,jj++)
        *ii= mod_diff(*ii,*jj);
    gb_fptr= &A[54];
    return A[55];
}

void gb_init_rand(long seed)
{
    long i;
    long prev= seed,next= 1;
    seed= prev= mod_diff(prev,0);
    A[55]= prev;
    for(i= 21;i;i= (i+21)%55){
        A[i]= next;
        next= mod_diff(prev,next);
        if(seed&1)seed= 0x40000000+(seed>>1);
        else seed>>= 1;
        next= mod_diff(next,seed);
        ;
        prev= A[i];
    }
    (void)gb_flip_cycle();
    (void)gb_flip_cycle();
    (void)gb_flip_cycle();
    (void)gb_flip_cycle();
    (void)gb_flip_cycle();
}

long gb_unif_rand(long m)
{
    unsigned long t= two_to_the_31-(two_to_the_31%m);
    long r;
    do{
        r= gb_next_rand();
    }while(t<=(unsigned long)r);
    return r%m;
}

//
// SHUFFLE ARRAY
//
void shuffle(int* array, int n)
{
    int i,j,t;
    for (j=n-1; j>0; j--) {
        i = gb_unif_rand(j);
        t = array[j];
        array[j] = array[i];
        array[i] = t;
    }
}

//
// HEAP
//
void minheapify(int *A, int *posinA, double *value, int pos, int numinA)
{
    int maxpos;
    int father, maxchild, rightchild;
    int fatherpos, maxchildpos, rightchildpos;
    maxpos = numinA - 1;
    fatherpos = pos;
    maxchildpos = 2*fatherpos + 1;
    maxchild = A[maxchildpos];
    while (maxchildpos<=maxpos) {
        father = A[fatherpos];
        if (maxchildpos<maxpos) {
            rightchildpos = maxchildpos + 1;
            rightchild = A[rightchildpos];
            if (value[maxchild]>value[rightchild]) {
                maxchild = rightchild;
                maxchildpos = rightchildpos;
            }
        }
        if (value[father]>value[maxchild]) {
            A[fatherpos] = maxchild;
            posinA[maxchild] = fatherpos;
            A[maxchildpos] = father;
            posinA[father] = maxchildpos;
            fatherpos = maxchildpos;
            maxchildpos = 2*fatherpos + 1;
            if (maxchildpos<=maxpos) maxchild = A[maxchildpos];
        }
        else maxchildpos = maxpos + 1;
    }
}

void decreasekey(int *A, int *posinA, double *value, int pos)
{
    int father, child;
    int fatherpos, childpos;
    if (pos>0) {
        childpos = pos;
        fatherpos = (childpos-1)/2;
        while (fatherpos>=0) {
            father = A[fatherpos];
            child = A[childpos];
            if (value[father]>value[child]) {
                A[fatherpos] = child;
                posinA[child] = fatherpos;
                A[childpos] = father;
                posinA[father] = childpos;
            }
            childpos = fatherpos;
            if (childpos>0) fatherpos = (childpos-1)/2;
            else fatherpos = -1;
        }
    }
}

//==================================================================================
//
// Heuristics for the linear ordering problem
// Note that the heuristics require a LOP in normal form, i.e. at most one of
// c(i,j) and c(j,i) is nonzero, and then it is positive.
// As a consequence the full matrix has to be stored and every ordering has
// a nonnegative value.
//
//==================================================================================
#define c(i,j) iomat[dv[i]+j]

int rgen1(int n);
int loCheck(int n,int CurVal,int *iomat,int *dv,int *ObjAt);
void updateBest(int n,int CurVal,int *ObjAt,int *PosOf,int *BestVal,int *BestObjAt,int *BestPosOf);
int loEval(int n,int *iomat,int *dv,int *ObjAt);
int randomOrder(int n,int *iomat,int *dv,int *ObjAt,int *PosOf,int *Val);
int locEnu(int n,int width,int *iomat,int *dv,int *ObjAt,int *PosOf,int *CurVal);
int linKer2(int n,int *iomat,int *dv,int *ObjAt,int *PosOf,int *CurVal);
int insertion(int n,int *iomat,int *dv,int *ObjAt,int *PosOf,int *CurVal);
void shuffleOrder(int n,int howmany,int *ObjAt,int *PosOf);
int lopheu();

int lopheu()
{
    if (dofixing || numcomponents > 1 || wnvar < 10 || nnod != numbottomnodes)
        return 0;
    int OscmVal;
    int NonZeros;
    
    int n;
    int t,h,a;
    //int doprint,dostat;
    int DiagSum, TotOrgSum; //TotSum;
    int i,j,k,min;
    
    if (pril) {
        printf("\nInitial heuristic ");fflush(stdout);
    }
    
    n = nnod;
    int BestVal = -10000000;
    int CurVal = 0;
    
    int ProbConst2 = objconst; // From oscm
    
    // Allocate
    
    int* ObjAt = (int *) malloc( (n+1) * sizeof(int) );
    int* PosOf = (int *) malloc( (n+1) * sizeof(int) );
    int* BestObjAt = (int *) malloc( (n+1) * sizeof(int) );
    int* BestPosOf = (int *) malloc( (n+1) * sizeof(int) );
    
    // Set dope vector and in initialize matrix
    
    int* iomat = (int *) malloc( (n*n+1) * sizeof(int) );
    int* dv = (int *) malloc( (n+1) * sizeof(int) );
    for (i=1;i<=n;i++) {
        dv[i] = (i-1)*n;
        for (j=1;j<=n;j++) c(i,j) = 0;
    }
    
    // Convert oscm LOP to maximization problem
    
    TotOrgSum = DiagSum = NonZeros = 0;
    for (a=0; a<wnvar; a++) {
        t = wtail[a] + 1;
        h = whead[a] + 1;
        k = (int) obj[a];
        if (k!=0) NonZeros++;
        k = -k;
        c(t,h) = k;
        TotOrgSum += k;
        if (t==h) DiagSum += k;
    }
    
    if (pril>1) printf("Number of sectors:  %10d\n",n);
    if (pril>1) printf("Original sum:       %10d\n",TotOrgSum);
    if (pril>1) printf("Diagonal sum:       %10d\n",DiagSum);
    if (pril>1) printf("Nonzeros:           %10d\n",NonZeros);
    
    // Transform matrix to normal form
    
    int AddConst = 0;
    for (i=1;i<n;i++) {
        c(i,i) = 0;
        for (j=i+1;j<=n;j++) {
            if (c(i,j)<c(j,i)) {
                min = c(i,j);
            }
            else {
                min = c(j,i);
            }
            AddConst += min;
            c(i,j) = c(i,j) - min;
            c(j,i) = c(j,i) - min;
        }
    }
    c(n,n) = 0;
    
    if (pril>1) printf("Additive constant:  %10d\n",AddConst);
    if (pril>1) printf("Oscm constant:      %10d\n",ProbConst2);
    
    //UpperBound = TotSum = TotOrgSum - DiagSum - 2*AddConst;
    //if (HPril) printf("\nTotal new sum:      %10d\n",TotSum);
    
    int Tries = 20;
    
    // LinKernighan + Shuffle
    
    if (pril>1) printf("\n");
    randomOrder(n,iomat,dv,ObjAt,PosOf,&CurVal);
    OscmVal = ProbConst2-(CurVal+AddConst);
    if (pril>0) printf("Random:             %10d [oscm %6d]\n",CurVal,OscmVal);
    
    loCheck(n,CurVal,iomat,dv,ObjAt);
    for (i=1;i<=Tries;i++) {
        // gb_unif_rand(100);
        if (pril) {
            printf("(%d) ",i);fflush(stdout);
        }
        linKer2(n,iomat,dv,ObjAt,PosOf,&CurVal);
        OscmVal = ProbConst2-(CurVal+AddConst);
        updateBest(n,CurVal,ObjAt,PosOf,&BestVal,BestObjAt,BestPosOf);
        if (pril>1) printf("Lin-Kernighan:      %10d [oscm %6d]\n",CurVal,OscmVal);
        loCheck(n,CurVal,iomat,dv,ObjAt);
        if (i<Tries) {
            shuffleOrder(n,(int) n/4, ObjAt, PosOf);
            CurVal = loEval(n,iomat,dv,ObjAt);
            OscmVal = ProbConst2-(CurVal+AddConst);
            if (pril>1) printf("Shuffle:            %10d [oscm %6d]\n",CurVal,OscmVal);
        }
    }
    
    OscmVal = ProbConst2-(BestVal+AddConst);
    
    double value = objconst;
    double bestvalue = objconst;
    for (int a=0; a<wnvar; a++) {
        int t = wtail[a];
        int h = whead[a];
        if (PosOf[t+1]<PosOf[h+1]) value += obj[a];
        if (BestPosOf[t+1]<BestPosOf[h+1]) bestvalue += obj[a];
    }
    if (pril>1) printf("\nCurrent value: %15d (%d)\n",(int) value,(int) (value-objconst));
    if (pril>1) printf("Best value:    %15d (%d)\n",(int) bestvalue,(int) (bestvalue-objconst));
    
    for (i=0; i<nnod; i++) {
        int boa = BestObjAt[i+1] - 1;
        activity[i] = boa;
        rank[boa] = i;
    }
    activityrank_num_crossings = bestvalue;
    if (pril>1) {
        printf("\ninitial ordering: ");
        for (i=0; i<nnod; i++) printf(" %d:%d",rank[activity[i]],activity[i]);
        printf("\n");
        printf("with %d crossings\n\n",(int) (activityrank_num_crossings+eps));
    }
    if (pril) printf("%d crossings\n",(int) (activityrank_num_crossings+eps));
    
    free(iomat);
    free(dv);
    free(ObjAt);
    free(PosOf);
    free(BestObjAt);
    free(BestPosOf);
    
    return 0;
}//lopheu

int loCheck(
    int n,
    int CurVal,
    int *iomat,
    int *dv,
    int *ObjAt
)
// Consistency check
{
    int *used;
    int i,j;

    used = (int *) malloc( (n+1) * sizeof(int) );
    for (i=1;i<=n;i++) used[i] = 0;
    for (i=1;i<=n;i++) used[ObjAt[i]] = 1;
    for (i=1;i<=n;i++) {
        if (used[i]==0) {
            for (j=1;j<=n;j++) printf("%d ",ObjAt[j]);
            printf("\n");
            printf("ERROR: Ordering incomplete! (%d missing)\n",i);
            exit(100);
        }
    }
    if(CurVal!=loEval(n,iomat,dv,ObjAt)) {
        printf("ERROR: Inconsistency of objective and ordering (%d!=%d)\n",
        CurVal,loEval(n,iomat,dv,ObjAt));
        exit(100);
    }
    free(used);
    return 0;
}//loCheck

void updateBest(
    int n,
    int CurVal,
    int *ObjAt,
    int *PosOf,
    int *BestVal,
    int *BestObjAt,
    int *BestPosOf
    )
{
    int i;

    if (CurVal>*BestVal) {
        *BestVal = CurVal;
        for (i=1;i<=n;i++) BestObjAt[i] = ObjAt[i];
        for (i=1;i<=n;i++) BestPosOf[i] = PosOf[i];
    }

}//updateBest

int loEval(
    int n,
    int *iomat,
    int *dv,
    int *ObjAt
    )
// Evaluate ordering given by ObjAt
{
    int i,j,sum;

    sum = 0;
    for (i=1;i<n;i++) {
        for (j=i+1;j<=n;j++) {
            sum += c(ObjAt[i],ObjAt[j]);
        }
    }
    return(sum);

}//loEval

int rgen1(
    int n
)
// Return random number between 1 and n
{
    return(gb_unif_rand(n) + 1);
}//rgen1

void shuffleOrder(
    int n,
    int howmany,
    int *ObjAt,
    int *PosOf
)
//Shuffle current order (a little bit)
{
    int i,t,o,p1,p2;
    for (t=1;t<=howmany;t++) {
        p1 = rgen1(n);
        p2 = rgen1(n);
        if (p1!=p2) {
                    o = ObjAt[p1];
                    ObjAt[p1] = ObjAt[p2];
                ObjAt[p2] = o;
                }
    }
        for (i=1;i<=n;i++) PosOf[ObjAt[i]] = i;
}//shuffleorder

int randomOrder(
    int n,
    int *iomat,
    int *dv,
    int *ObjAt,
    int *PosOf,
    int *Val
    )
// Generate a random ordering
{
    int i,j,k,t;
    int sum,totsum;
    for (i=1;i<=n;i++) ObjAt[i] = i;
    for (i=1;i<=n;i++) {
        j = rgen1(n);
        k = rgen1(n);
        t = ObjAt[j];
        ObjAt[j] = ObjAt[k];
        ObjAt[k] = t;
    }
    for (i=1;i<=n;i++) PosOf[ObjAt[i]] = i;

    totsum = sum = 0;
    for (i=1;i<n;i++) {
        for (j=i+1;j<=n;j++) {
            sum += c(ObjAt[i],ObjAt[j]);
            totsum += c(ObjAt[i],ObjAt[j]) + c(ObjAt[j],ObjAt[i]);
        }
    }

    if (sum<totsum/2+1) {
        for (i=1;i<=n/2;i++) {
            t = ObjAt[i];
            ObjAt[i] = ObjAt[n-i+1];
            ObjAt[n-i+1] = t;
        }
        for (i=1;i<=n;i++) PosOf[ObjAt[i]] = i;
        *Val = totsum - sum;

    }
    else {
        *Val = sum;
    }
    return 0;
}//randomOrder

int linKer2(
    int n,
    int *iomat,
    int *dv,
    int *ObjAt,
    int *PosOf,
    int *CurVal
    )
// Perform Kernighan-Lin type heuristic (2)
{
    int progress;
    int value,gain,maxgain,totalgain,bestgain,bestinsgain;
    int i,j,l,nimove,v,w;
    int newj,bestjpos,bestinsj,bestv,bestnewj,maxpos;
    int *actold,*available,*movej,*movenewj,*movev;
    bestv = -1;
    bestinsj = -1;
    bestnewj = -1;
    bestjpos = -1;
    maxpos = -1;

    actold = (int *) malloc( (n+1) * sizeof(int) );
    available = (int *) malloc( (n+1) * sizeof(int) );
    movej = (int *) malloc( (n+1) * sizeof(int) );
    movenewj = (int *) malloc( (n+1) * sizeof(int) );
    movev = (int *) malloc( (n+1) * sizeof(int) );

    value = *CurVal;

    do {
        progress = 0;
        for (i=1;i<=n;i++) {
            available[i] = 1;
            actold[i] = ObjAt[i];
        }

        maxgain = -value;
        totalgain = 0;

        for (nimove=1;nimove<=2*n/3;nimove++) {
            bestinsgain = -value;
            for (j=1;j<=n;j++) if (available[ObjAt[j]]) {
                v = ObjAt[j];

                // Find best position for v (change position in any case)
                // a) positions before j
                bestgain = -value;
                for (newj=1;newj<j;newj++) {
                    gain = 0;
                    for(i=newj;i<j;i++) {
                        w = ObjAt[i];
                        gain += c(v,w) - c(w,v);
                    }
                    if (gain>bestgain) {
                        bestjpos = newj;
                        bestgain = gain;
                    }
                }
                // b) positions after j
                for (newj=j+1;newj<=n;newj++) {
                    gain = 0;
                    for(i=j+1;i<=newj;i++) {
                        w = ObjAt[i];
                        gain += c(w,v) - c(v,w);
                    }
                    if (gain>bestgain) {
                        bestjpos = newj;
                        bestgain = gain;
                    }
                }

                // Store best possiblity over all j
                if (bestgain>bestinsgain) {
                    bestinsj = j;
                    bestv = v;
                    bestnewj = bestjpos;
                    bestinsgain = bestgain;
                }
            }

            if (bestv == -1 || bestinsj == -1 || bestnewj == -1) {
                if (pril) printf("bestv/bestinsj/bestnewj not set\n");
                return 0;
            }

            // Tentatively realize best insertion
            available[bestv] = 0;
            movej[nimove] = bestinsj;
            movenewj[nimove] = bestnewj;
            movev[nimove] = bestv;

            if (bestnewj<bestinsj) {
                for (l=bestinsj;l>bestnewj;l--)
                    ObjAt[l] = ObjAt[l-1];
                ObjAt[bestnewj] = bestv;
            }
            else {
                for (l=bestinsj;l<bestnewj;l++)
                    ObjAt[l] = ObjAt[l+1];
                ObjAt[bestnewj] = bestv;
            }

            totalgain += bestinsgain;
            if (totalgain>maxgain) {
                maxgain = totalgain;
                maxpos = nimove;
            }
        }

        // First undo all moves
        for (i=1;i<=n;i++) ObjAt[i] = actold[i];
        
        if (maxgain>0) {
            // Permanently realize moves up to maxpos
            progress = 1;
            value += maxgain;
            if (maxpos == -1) {
                if (pril) printf("maxpos not set\n");
                return 0;
            }
            for (nimove=1;nimove<=maxpos;nimove++) {
                if (movenewj[nimove]<movej[nimove]) {
                    for (l=movej[nimove];l>movenewj[nimove];l--)
                        ObjAt[l] = ObjAt[l-1];
                    ObjAt[movenewj[nimove]] = movev[nimove];
                }
                else {
                    for (l=movej[nimove];l<movenewj[nimove];l++)
                        ObjAt[l] = ObjAt[l+1];
                    ObjAt[movenewj[nimove]] = movev[nimove];
                }
            }
        }
    } while (progress);

    for (i=1;i<=n;i++) {
        PosOf[ObjAt[i]] = i;
    }
    *CurVal = value;

    free(actold);
    free(available);
    free(movej);
    free(movenewj);
    free(movev);
    return 0;
}//linKer2
#undef c

//
// EXPLOITATION OF LP SOLUTIONS
//
int asexploitlp(const double* wxval, int cbcCallback, double& cbcSolutionValue, double* betterSolution)
{
    double starttime = processtime();

    int i,j,v,w,a;
    int numinPQ;
    int mininvalnode;
    int neighbor, neighborarc;
        
    double  locvalue;
    
    double *inval;
    int *indeg;
    
    int   *locact;
    int   *locran;
    int   *PQ;
    int   *posinPQ;
    int   *stack;
    
    locact  = (int *) malloc(nnod*sizeof(int));
    locran  = (int *) malloc(numbottomnodes*sizeof(int));

    // count fractional arcs
    int isintegral = true;
    for (int i = 0; i < wnvar; i++) {
        if (wxval[i] > eps && wxval[i] < omeps) {
            isintegral = false;
            break;
        }
    }

    // first try exact topsort
    int solutionfound = false;
    if (isintegral) {
        indeg = (int *) malloc(numbottomnodes*sizeof(int));

        for (i=0; i<nnod; i++) indeg[node[i]] = 0;
        for (a=0; a<wnvar; a++) {
            if (wxval[a] >= omeps)
                indeg[whead[a]]++;
            else
                indeg[wtail[a]]++;
        }
        for (a=0; a<numfixed; a++) {
            indeg[fhead[a]]++;
        }

        stack = (int*) malloc(numbottomnodes*sizeof(int));
        int stacksize = 0;
        // add zero indeg nodes to stack
        for (i=0; i<nnod; i++) {
            if (indeg[node[i]] == 0) {
                stack[stacksize++] = node[i];
            }
        }
        int nodecount = 0;
        while (stacksize > 0) {
            v = stack[stacksize - 1];
            stacksize--;

            locran[v] = nodecount++;

            // adjust indeg of neighbors for outgoing warcs (wxval)
            neighborarc = wfirst[v];
            while (neighborarc!=-1) {
                neighbor = whead[neighborarc];
                if (wxval[neighborarc] >= omeps) {
                    indeg[neighbor]--;
                    if (indeg[neighbor] == 0)
                        stack[stacksize++] = neighbor;
                }
                neighborarc = wnext[neighborarc];
            }

            // adjust indeg of neighbors for ingoing warcs (1-wxval)
            neighborarc = wfirstin[v];
            while (neighborarc!=-1) {
                neighbor = wtail[neighborarc];
                if (wxval[neighborarc] < omeps) {
                    indeg[neighbor]--;
                    if (indeg[neighbor] == 0)
                        stack[stacksize++] = neighbor;
                }
                neighborarc = wnextin[neighborarc];
            }
            
            // adjust indeg of neighbors for fixed arcs
            neighborarc = ffirst[v];
            while (neighborarc!=-1) {
                neighbor = fhead[neighborarc];
                indeg[neighbor]--;
                if (indeg[neighbor] == 0)
                    stack[stacksize++] = neighbor;
                neighborarc = fnext[neighborarc];
            }

            if (nodecount > nnod) {
                printf("ERROR: nodecount>nnod\n");
                exit(1);
            }
        }
        // check whether topsort is complete
        if (nodecount == nnod) {
            solutionfound = true;
            if (pril) printf("ASPEXPLOITLP Solution found using exact topsort.\n");
        }

        free(indeg);
        free(stack);
    }

    // heuristic topsort
    if (!solutionfound) {
        inval   = (double *) malloc(numbottomnodes*sizeof(double));
        PQ      = (int *) malloc(nnod*sizeof(int));
        posinPQ = (int *) malloc(numbottomnodes*sizeof(int));

        for (i=0; i<nnod; i++) inval[node[i]] = 0.0;
        for (a=0; a<wnvar; a++) {
            inval[whead[a]] += wxval[a];
            inval[wtail[a]] += 1 - wxval[a];
        }
        for (a=0; a<numfixed; a++) {
            inval[fhead[a]] += 2.0 * wnvar; // high weight to fixed edges
        }
        
        // make a priority queue according to inval
        for (i=0; i<nnod; i++) {
            PQ[i] = node[i];
            posinPQ[node[i]] = i;
        }
        numinPQ = nnod;
        for (i=numinPQ/2-1; i>=0; i--) {
            minheapify(PQ,posinPQ,inval,i,numinPQ);
        }
            
        // rank all nodes
        for (i=0; i<nnod; i++) {
            mininvalnode = PQ[0];
            locran[mininvalnode] = i;
            // restore heap
            PQ[0] = PQ[--numinPQ];
            posinPQ[PQ[0]] = 0;
            if (numinPQ>1) {
                minheapify(PQ,posinPQ,inval,0,numinPQ);
            }
            
            // adjust invalue of neighbors for outgoing warcs (wxval)
            neighborarc = wfirst[mininvalnode];
            while (neighborarc!=-1) {
                neighbor = whead[neighborarc];
                inval[neighbor] -= wxval[neighborarc];
                if ((numinPQ>1)&&(posinPQ[neighbor]<numinPQ)) {
                    decreasekey(PQ,posinPQ,inval,posinPQ[neighbor]);
                }
                neighborarc = wnext[neighborarc];
            }

            // adjust invalue of neighbors for ingoing warcs (1-wxval)
            neighborarc = wfirstin[mininvalnode];
            while (neighborarc!=-1) {
                neighbor = wtail[neighborarc];
                inval[neighbor] -= 1 - wxval[neighborarc];
                if ((numinPQ>1)&&(posinPQ[neighbor]<numinPQ)) {
                    decreasekey(PQ,posinPQ,inval,posinPQ[neighbor]);
                }
                neighborarc = wnextin[neighborarc];
            }
            
            // adjust invalue of neighbors for fixed arcs
            neighborarc = ffirst[mininvalnode];
            while (neighborarc!=-1) {
                neighbor = fhead[neighborarc];
                inval[neighbor] -= 2.0 * wnvar;
                if ((numinPQ>1)&&(posinPQ[neighbor]<numinPQ)) {
                    decreasekey(PQ,posinPQ,inval,posinPQ[neighbor]);
                }
                neighborarc = fnext[neighborarc];
            }
        }
        free(inval);
        free(PQ);
        free(posinPQ);
    }

    for (i=0; i<nnod; i++) locact[locran[node[i]]] = node[i];

    // compute value
    locvalue = 0.0;
    for (i=0; i<nnod; i++) {
        v = locact[i];
        for (j=i+1; j<nnod; j++) {
            w = locact[j];
            locvalue += ((double) compcrosstable[poscct[v]*nnod+poscct[w]]);
        }
    }

    // better solution found than stored in activity/rank? -> update
    if (locvalue < activityrank_num_crossings) {
        for (i=0; i<nnod; i++) activity[i] = locact[i];
        for (i=0; i<numbottomnodes; i++) rank[i] = locran[i];
        activityrank_num_crossings = locvalue;
        if (!cbcCallback && pril>1) {
            printf("ASPEXPLOITLP upper bound improved to %.0lf crossings.\n",locvalue);
        }
    }

    // give better solution to cbc
    int improvedcbcsolution = 0;
    if (cbcCallback) {
        // compute obj value
        double objval = 0.0;
        for (a = 0; a < wnvar; a++) {
            int head = whead[a];
            int tail = wtail[a];
            if (locran[tail] < locran[head])
                objval += obj[a];
        }
        if (objval + objconst != locvalue) {
            printf("Error: objval + objconst != locvalue %f %f\n", objval + objconst, locvalue);
            exit(1);
        }

        if (objval < cbcSolutionValue) {
            // store new solution in cbcSolutionValue and betterSolution
            cbcSolutionValue = objval;
            for (a = 0; a < wnvar; a++) {
                int head = whead[a];
                int tail = wtail[a];
                if (locran[tail] < locran[head])
                    betterSolution[a] = 1.0;
                else
                    betterSolution[a] = 0.0;
            }
            if (pril) printf("ASPEXPLOITLP CBC upper bound improved to %.0lf crossings. (obj: %.0f)\n",locvalue, objval);

            improvedcbcsolution = 1;
        }
    }
    
    free (locact);
    free (locran);

    exploitationtime += processtime() - starttime;
    
    return improvedcbcsolution;
}


//
// CYCLE SEPARATION
//
int fascycle(Graph* g, int onlyfractional, OsiCuts& osiCuts)
{
    if (pril) {
        printf("FASCYCLE"); fflush(stdout);
    }
    
    /* local scalars */
    int ngen;
    int v;
    int e;
    int i;
    int t,h;
    int minv = -1;
    int last;
    int numcyclearcs;
    int numcyclearcs_total;
    int numinPQ;
    int neighborarc;
    int checknode,checknodepos;
    int predarcchecknode;
    int posh;
    int pq0;

    double infty;
    double distminv;
    double distchecknode;
    double lengthneighborarc;
    double cyclelength;
    double w;
    
    int *PQ;
    int *posinPQ;
    int *pred;
    int *predarc;
    double *dist;

    double rhs;

    double eps = 0.0001;      // epsilon
    double omeps = 1.0 - eps; // 1 - epsilon
    
    // length = (double*) malloc(nvar*sizeof(double));
    dist = (double*) malloc(numbottomnodes*sizeof(double));
    PQ = (int*) malloc(numbottomnodes*sizeof(int));
    posinPQ = (int*) malloc(numbottomnodes*sizeof(int));
    pred = (int*) malloc(numbottomnodes*sizeof(int));
    predarc = (int*) malloc(numbottomnodes*sizeof(int));

    /*---------- body of fascycle ---------*/
    
    /* initialize */
    ngen = 0;
    infty = 1e10;
    //double edgeEps = fmin(1 / (nedge + 1), eps) / 100;
    // int* edgePermutation = (int*)malloc(sizeof(int) * nedge);
    // for (int i = 0; i < nedge; i++)
    //     edgePermutation[i] = i;
    // shuffle(edgePermutation, nedge);
    for (e=0; e < g->nedges; e++) {
        // g->length[e] = fmin(fmax(1.0 - g->xval[e], 0.0), 1.0);// + edgeEps;
        g->unused[e] = 1;
    }
    /* look for violated cycles */
    for (e = 0; e < g->nedges; e++) if (g->unused[e]) {
        // int firstarc = edgePermutation[e];
        if (onlyfractional && (g->fixed[e] || g->xval[e] > omeps || g->xval[e] < eps))
            continue;
        // printf("%f\n", g->length[e]);

        // printf("\nChecking edge %d->%d\n", g->tail[e], g->head[e]);
        t = g->tail[e];
        h = g->head[e];
        // printf("Edge %d->%d (%d)\n", t, h, e);

        w = g->fixed[e] ? 0 : (1 - g->xval[e]);
        if (w<omeps) { // might consider smaller as heuristic measure
            // g->length[e] = infty;
            /* find shortest h-t-path or conclude that there is no h-t-path (Dijkstra) */
            for (i=0; i<nnod; i++) {
                v = node[i];
                dist[v] = infty;
                pred[v] = -1;
            }
            dist[h] = 0.0;
            /* initialize priority queue */
            for (i=0; i<nnod; i++) {
                v = node[i];
                PQ[i] = v;
                posinPQ[v] = i;
            }
            posh = posinPQ[h];
            pq0 = PQ[0];
            PQ[posh] = PQ[0];
            PQ[0] = h;
            posinPQ[h] = 0;
            posinPQ[pq0] = posh;
            numinPQ = nnod;
            /* Dijkstra */
            while (numinPQ) {
                /* extractmin */
                minv = PQ[0];
                distminv = dist[minv];
                /* restore priority queue */
                numinPQ--;
                last = PQ[numinPQ];
                PQ[0] = last;
                posinPQ[last] = 0;
                PQ[numinPQ] = minv;
                posinPQ[minv] = numinPQ;
                if (numinPQ>1) minheapify(PQ,posinPQ,dist,0,numinPQ);
                /* adjust priorities of neighbors */
                if (numinPQ) {
                    neighborarc = g->firstout[minv];
                    while (neighborarc!=-1) {
                        // printf("%d\n", neighborarc);
                        checknode = g->head[neighborarc];
                        // if (t == 27 && h == 254 && checknode == 27) {
                        //     printf("checknode %d\n", checknode);
                        // }
                        checknodepos = posinPQ[checknode];
                        if ((checknodepos>=0)&&(checknodepos<numinPQ)&&neighborarc != w) {
                            distchecknode = dist[checknode];
                            lengthneighborarc = g->fixed[neighborarc] ? 0 : (1 - g->xval[neighborarc]);
                            if (distminv+lengthneighborarc<distchecknode) {
                                // if (t == 27 && h == 254) {
                                //     printf("%d -> %d (%.8f to %.8f)\n", minv, checknode, distchecknode, distminv+lengthneighborarc);
                                // }
                                distchecknode = distminv+lengthneighborarc;
                                dist[checknode] = distchecknode;
                                pred[checknode] = minv;
                                predarc[checknode] = neighborarc;

                                /* adjust priority queue */
                                decreasekey(PQ,posinPQ,dist,checknodepos);
                            }
                        } /* if (checknodepos>=0)&&(checknodepos<numinPQ) */
                        neighborarc = g->nextout[neighborarc];
                    } /* while (neighborarc!=-1) */
                } /* if (numinPQ) */
                if ((dist[PQ[0]]>=infty)||(minv==t)) numinPQ = 0;
            } /* while (numinPQ) */
            /* restore length of e */
            // g->length[e] = w;

            if (minv==t) {
                /* compute length of cycle */
                numcyclearcs_total = 0;
                numcyclearcs = 0;
                cyclelength = 0.0;
                rhs = 0.0;

                if (!g->fixed[e]) {
                    // not fixed, insert edge to constraint
                    cind[numcyclearcs] = g->warc[e];
                    if (g->tail[e] < g->head[e]) {
                        cval[numcyclearcs] = 1.0;
                    } else {
                        cval[numcyclearcs] = -1.0;
                        rhs -= 1.0;
                    }
                    cyclelength += w;
                    numcyclearcs++;
                }
                graph_cind[numcyclearcs_total++] = e;

                checknode = t;
                while (checknode!=h) {
                    predarcchecknode = predarc[checknode];

                    if (!g->fixed[predarcchecknode]) {
                        // not fixed, insert edge to constraint
                        cind[numcyclearcs] = g->warc[predarcchecknode];
                        if (g->tail[predarcchecknode] < g->head[predarcchecknode]) {
                            cval[numcyclearcs] = 1.0;
                        } else {
                            cval[numcyclearcs] = -1.0;
                            rhs -= 1.0;
                        }
                        if (!g->fixed[predarcchecknode])
                            cyclelength += 1 - g->xval[predarcchecknode];
                        numcyclearcs++;
                    }
                    graph_cind[numcyclearcs_total++] = predarcchecknode;

                    checknode = pred[checknode];
                }

                /* generate violated inequality */
                if (cyclelength < omeps && numcyclearcs_total > 3) {
                    rhs += numcyclearcs - 1;
                    // printf("adding constraint with rhs %.2f\n",rhs);
                    OsiRowCut rc;
                    rc.setLb(-COIN_DBL_MAX);
                    rc.setUb(rhs);
                    rc.setRow(numcyclearcs, cind, cval);
                    rc.setGloballyValid(true);
                    osiCuts.insert(rc);

                    ngen++;
                    
                    for (i=0; i<numcyclearcs_total; i++) g->unused[graph_cind[i]] = 0;
                    if (ngen >= max_gen_cycles) {
                        goto enough;
                    }
                } /* if (cyclelength < omeps) */
            } /* if (minv==t) */
        } /* if (w<eps) */
        g->unused[e] = 0;
    } /* for (e=0; e<nvar; e++) */

enough:
    if (pril > 0) printf(" %d violated cycles generated.\n",ngen);

    free(dist);
    free(PQ);
    free(posinPQ);
    free(pred);
    free(predarc);

    numfraccyclecuts += ngen;

    return ngen;
} /* end fascycle */

//
// MOEBIUS SEPARATION
//
int moebiusSep(const double* wxval, OsiCuts& osiCuts) {
    double startTime = processtime();
    // const double MAX_MOEBIUS_TIME = 200000;

    if (pril) {
        printf("MOEBIUS "); fflush(stdout);
    }

#define translate(i,j)                      \
    ori = origin[(i)*numbottomnodes+(j)];   \
    if (!ori) {                             \
        rhs -= xmat[i*numbottomnodes+j];    \
    }                                       \
    else if (ori>0.0) {                     \
        cind[cnvar] = ori - 1;              \
        cval[cnvar++] = 1.0;                \
        lhs += wxval[ori-1];                \
    }                                       \
    else {                                  \
        cind[cnvar] = -ori - 1;             \
        cval[cnvar++] = -1.0;               \
        rhs -= 1.0;                         \
        lhs -= wxval[-ori-1];               \
    }                                       \

    struct constr {
        int ccnvar;
        int ccind[11];
        double ccval[11];
        double crhs;
        double cviolation;
    };
    
    
    Graph F;
    int numfracarcs;
    int i,j,k,v,w,a,t,h;
    int mb,mt,lb,lt,rb,rt;
    int v1,v2,mu,a01,a12;
    int gs1,gs2;
    int ge1,ge2;
    int g1,g2;
    int p1,p2;
    int numpaths;
    int numgroups;
    long long numcandidates;
    const long long MAX_NUM_CANDIDATES = 5000000000; // 5e9
    int cnvar;
    int ori;
    int ind;
    int ngen = 0;
    int numstored = 0;
    int maxnumstored = 600000;

    double x;
    double f0,f1,f2,fv1mb, fv2mt;
    double rhs;
    double val;
    double vala,valb;
    double rblt, lbrt;

    double *xmat;
    double lhs;
    
    int *groupstart;
    int *origin;
    int *cind;
    double *cval;

    int *av1, *av2;
    double *pathweight;
    int *randomarc;

    // don't run on very large instances
    if (numbottomnodes > 5000) {
        return 0;
    }

    // count fractional arcs
    numfracarcs = 0;
    for (a=0; a<wnvar; a++) {
        x = wxval[a];
        if ((x>eps)&&(x<omeps)) numfracarcs++;
    }
    numfracarcs *= 2;

    if (numfracarcs < 5) {
        return 0;
    }

    // allocations
    xmat = (double *) malloc(numbottomnodes*numbottomnodes*sizeof(double));
    origin = (int *) malloc(numbottomnodes*numbottomnodes*sizeof(int));
    av1 = (int *) malloc(numbottomnodes*numbottomnodes*sizeof(int));
    av2 = (int *) malloc(numbottomnodes*numbottomnodes*sizeof(int));
    pathweight = (double *) malloc(numbottomnodes*numbottomnodes*sizeof(double));
    groupstart = (int *) malloc((numbottomnodes+1)*sizeof(int));
    cind = (int *) malloc(11*sizeof(int));
    cval = (double *) malloc(11*sizeof(double));
    constr* constraints = (constr*)malloc(maxnumstored * sizeof(constr));

    // initialize xmat
    for (i=0; i<numbottomnodes*numbottomnodes; i++) {
        xmat[i] = -1.0;
    }
    
    randomarc = (int *) malloc(numfracarcs*sizeof(int));

    // make fractional graph and copy w-part of xmat
    allocate_graph(F,numbottomnodes,numfracarcs, 0);
    init_graph(F,numbottomnodes);
    for (a=0; a<wnvar; a++) {
        x = fmin(fmax(wxval[a], 0), 1);
        t = wtail[a];
        h = whead[a];
        if ((x>eps)&&(x<omeps)) {
            insert_wedge(F,t,h,a,x);
            insert_wedge(F,h,t,a,1.0-x);
        }
        xmat[t*numbottomnodes+h] = x;
        origin[t*numbottomnodes+h] = a+1;
        xmat[h*numbottomnodes+t] = 1.0-x;
        origin[h*numbottomnodes+t] = -a-1;
    }
    
    // copy fixed values to xmat
    for (a=0; a<numfixed; a++) {
        t = ftail[a];
        h = fhead[a];
        xmat[t*numbottomnodes+h] = 1.0;
        origin[t*numbottomnodes+h] = 0;
        xmat[h*numbottomnodes+t] = 0.0;
        origin[h*numbottomnodes+t] = 0;
    }
    
    // show xmat
    if (pril>1) {
        int numzeroes = 0;
        int numones = 0;
        int numfracs = 0;
        int numminusones = 0;
        for (v=0; v<numbottomnodes; v++) {
            for (w=0; w<numbottomnodes; w++) {
                x = xmat[v*numbottomnodes+w];
                if ((x>eps)&&(x<omeps)) numfracs++;
                else if (x<-0.9) numminusones++;
                else if (x>=omeps) numones++;
                else numzeroes++;
            }
        }
        printf("%d ones, %d zeros, %d fracs, %d minusones\n",numones,numzeroes,numfracs,numminusones);
    }
    
    // make random permutation of fractional arcs
    for (i=0; i<numfracarcs; i++) randomarc[i] = i;
    for (j=numfracarcs-1; j>0; j--) {
        i = gb_unif_rand(j);
        t = randomarc[j];
        randomarc[j] = randomarc[i];
        randomarc[i] = t;
    }
    
    // list promising configurations
    numcandidates = 0;
    for (k=0; k<numfracarcs; k++) {
        mu = randomarc[k];
        numpaths= 0;
        numgroups= 0;
        mb = F.tail[mu];
        mt = F.head[mu];
        f0 = F.xval[mu];
        a01 = F.firstout[mt];
        while (a01>=0) {
            f1 = F.xval[a01];
            v1 = F.head[a01];
            fv1mb = xmat[v1*numbottomnodes+mb];
            if (v1 != mb && v1 != mt && fv1mb != -1.0 /*arbitrary ordering*/) {
                groupstart[numgroups++] = numpaths;
                a12 = F.firstout[v1];
                while (a12>=0) {
                    f2 = F.xval[a12];
                    v2 = F.head[a12];
                    fv2mt = xmat[v2*numbottomnodes+mt];
                    if (v2 != mb && v2 != mt && fv2mt != -1.0 /*arbitrary ordering*/) {
                        valb = f0+f1+f2+fv1mb+fv2mt;
                        vala = 5.0 - valb;
                        if ((vala>=3.0)||(valb>=3.0)) {
                            av1[numpaths] = v1;
                            av2[numpaths] = v2;
                            pathweight[numpaths] = vala - (1-f0);
                            numpaths++;
                        }
                    } // if (usable[v2])
                    a12 = F.nextout[a12];
                } // while (a12>=0)
            } // if (usable[v1])
            a01 = F.nextout[a01];
        } // while (a01>=0)

        groupstart[numgroups] = numpaths;

        // examine suitable pairs
        for (g1=0; g1<numgroups-1; g1++) {
            gs1 = groupstart[g1];
            ge1 = groupstart[g1+1];
            for (g2=g1+1; g2<numgroups; g2++) {
                gs2 = groupstart[g2];
                ge2 = groupstart[g2+1];
                for (p1=gs1; p1<ge1; p1++) {
                    lb = av1[p1];
                    lt = av2[p1];
                    for (p2=gs2; p2<ge2; p2++) {
                        rb = av1[p2];
                        rt = av2[p2];
                        assert(lb != rb);
                        
                        rblt = xmat[rb*numbottomnodes+lt];
                        lbrt = xmat[lb*numbottomnodes+rt];
                        if (lt!=rt && lt != rb && rt != lb && rblt != -1.0 && lbrt != -1.0 /*arbitrary ordering*/) {

                            numcandidates++;
                            
                            if (pathweight[p1] + pathweight[p2] + (1-f0) +
                                rblt + lbrt > 8.0+eps) {
                            
                                cnvar = 0;
                                rhs = 8.0;
                                lhs = 0.0;
                                translate(mt,mb);
                                translate(lb,mt);
                                translate(lt,lb);
                                translate(rb,mt);
                                translate(rt,rb);
                                translate(mb,lb);
                                translate(mb,rb);
                                translate(mt,lt);
                                translate(rb,lt);
                                translate(mt,rt);
                                translate(lb,rt);
                                
                                assert(lhs > rhs);
                                
                                // store constraint
                                constraints[numstored].crhs = rhs;
                                constraints[numstored].cviolation = lhs - rhs;
                                constraints[numstored].ccnvar = cnvar;
                                for (i=0; i<cnvar; i++) {
                                    constraints[numstored].ccind[i] = cind[i];
                                    constraints[numstored].ccval[i] = cval[i];
                                }
                                numstored++;
                                if (numstored>=maxnumstored) {
                                    goto enough;
                                }
                            }
                            
                            numcandidates++;
                            
                            if (10.0 - pathweight[p1] - pathweight[p2] + f0
                                - rblt - lbrt > 8.0+eps) {
                                
                                cnvar = 0;
                                rhs = 8.0;
                                lhs = 0.0;
                                translate(mb,mt);
                                translate(mt,lb);
                                translate(lb,lt);
                                translate(mt,rb);
                                translate(rb,rt);
                                translate(lb,mb);
                                translate(rb,mb);
                                translate(lt,mt);
                                translate(lt,rb);
                                translate(rt,mt);
                                translate(rt,lb);

                                assert(lhs > rhs);
                            
                                // store constraint
                                constraints[numstored].crhs = rhs;
                                constraints[numstored].cviolation = lhs - rhs;
                                constraints[numstored].ccnvar = cnvar;
                                for (i=0; i<cnvar; i++) {
                                    constraints[numstored].ccind[i] = cind[i];
                                    constraints[numstored].ccval[i] = cval[i];
                                }
                                numstored++;
                                if (numstored>=maxnumstored) {
                                    goto enough;
                                }
                            }
                            
                            if (numcandidates>=MAX_NUM_CANDIDATES) {
                                goto enough;
                            }
                        }
                    }
                }
            }
        }
    } // for (k=0; k<numfracarcs; k++)
    
    enough:
    
    // shuffle constraints
    for (j=numstored-1; j>0; j--) {
        i = gb_unif_rand(j);
        constr t = constraints[j];
        constraints[j] = constraints[i];
        constraints[i] = t;
    }

    ngen = 0;
    int n = numstored; if (max_gen_cycles < n) n = max_gen_cycles;
    for (i = 0; i < n; i++) {
        // check inequality
        for (j=0; j<constraints[i].ccnvar; j++) {
            ind = constraints[i].ccind[j];
            val = constraints[i].ccval[j];
            if (ind < 0 || ind >= wnvar || (val != -1.0 && val != 1.0)) {
                printf("ERROR %d\n", ind);
                exit(1);
            }
        }
        
        // insert row
        OsiRowCut rc;
        rc.setLb(-COIN_DBL_MAX);
        rc.setUb(constraints[i].crhs);
        rc.setRow(constraints[i].ccnvar, constraints[i].ccind, constraints[i].ccval);
        rc.setGloballyValid(true);
        osiCuts.insert(rc);
        ngen++;
    }
    
    if (pril) printf(" %lld candidates, %d constraints to choose from, %d inequalities generated.\n",numcandidates,numstored,ngen);

    free(constraints);
    free_graph(F);
    free(xmat);
    free(origin);
    free(av1);
    free(av2);
    free(pathweight);
    free(groupstart);
    free(cind);
    free(cval);
    free(randomarc);

    separationtimemoebius += (processtime() - startTime);

    nummoebiuscuts += ngen;

    return ngen;
#undef translate
}

//
// HEURISTIC SEPARATION OF 3-CYCLES
// Note: DOESNT SUPPORT FIXED EDGES
//

int fast3CycleSep(Graph* g, OsiCuts& osiCuts) {
    if (pril) {
        printf("3-CYCLES   ");fflush(stdout);
    }
    // DFS to find B-edges

    int i, j, arc;
    int startNode, node;
    double rhs, cycleval;
    int b_edge_index;
    int b_arc, b_tail, b_head;
    int neighbor, current_neighborarc;

    // arrays for dfs
    int* visited = (int*)malloc(sizeof(int) * numbottomnodes);
    int* nodepermutation = (int*)malloc(sizeof(int) * numbottomnodes);
  
    // stacks for dfs
    int* node_stack = (int*) malloc(sizeof(int) * numbottomnodes);
    int* current_neighborarc_stack = (int*) malloc(sizeof(int) * numbottomnodes);
    int stack_size = 0;

    int* is_b_edge = (int*) malloc(sizeof(int) * g->nedges);
    int* b_edges = (int*) malloc(sizeof(int) * g->nedges);
    int n_b_edges = 0;

    for (i = 0; i < numbottomnodes; i++) {
        visited[i] = 0;
        nodepermutation[i] = i;
    }
    for (i = 0; i < g->nedges; i++) {
        is_b_edge[i] = 0;
    }

    shuffle(nodepermutation, numbottomnodes); // random node permutation
    for (i = 0; i < numbottomnodes; i++) { // restart dfs on each node
        startNode = nodepermutation[i];
        if (visited[startNode] != 0)
            continue;
        node_stack[stack_size] = startNode; // add start node to stack
        current_neighborarc_stack[stack_size++] = g->firstout[startNode];

        while (stack_size > 0) {
            node = node_stack[stack_size - 1];
            if (visited[node] == 0) { // newly discovered
                visited[node] = 1;
            }
            // traverse all neighbors
            current_neighborarc = current_neighborarc_stack[stack_size - 1];
            if (current_neighborarc != -1) {
                current_neighborarc_stack[stack_size - 1] = g->nextout[current_neighborarc];
                neighbor = g->head[current_neighborarc];
                if (visited[neighbor] == 1) {
                    // cycle found
                    b_edges[n_b_edges] = current_neighborarc;
                    is_b_edge[current_neighborarc] = 1;
                    n_b_edges++;
                } else if (visited[neighbor] == 0) {
                    // add to stack
                    node_stack[stack_size] = neighbor;
                    current_neighborarc_stack[stack_size++] = g->firstout[neighbor];
                }
            } else {
                visited[node] = 2;
                stack_size--;
            }
            if (stack_size > numbottomnodes) {
                printf("stack size %d>=nnod (%d)\n", stack_size, numbottomnodes);
                exit(1);
            }
        }
    }

    // shuffle b-edges
    shuffle(b_edges, n_b_edges);

    int* in_edge = (int*) malloc(sizeof(int) * numbottomnodes);

    struct constr {
        int cind[3];
        double cval[3];
        double rhs;
        double violation;
    };
    constr t;
    int MAX_STORED_CONSTRAINTS = 600000;
    constr* constraints = (constr*)malloc(MAX_STORED_CONSTRAINTS * sizeof(constr));
    int nconstr = 0;

    for (b_edge_index = 0; b_edge_index < n_b_edges; b_edge_index++) {
        for (i = 0; i < numbottomnodes; i++)
            in_edge[i] = -1;
        b_arc = b_edges[b_edge_index];
        b_tail = g->tail[b_arc];
        b_head = g->head[b_arc];

        // store edges (b_head, u) in in_edge[u]
        arc = g->firstout[b_head];
        while (arc != -1) {
            if (!is_b_edge[arc]) // to prevent duplicate 3-cycles
                in_edge[g->head[arc]] = arc;
            arc = g->nextout[arc];
        }
        // compare with incoming edges of b_tail
        arc = g->firstin[b_tail];
        while (arc != -1) {
            neighbor = g->tail[arc];
            if (in_edge[neighbor] != -1) {
                // 3-cycle found
                rhs = 2;
                cycleval = 0.0;

                cind[0] = g->warc[b_arc];
                if (g->tail[b_arc] < g->head[b_arc]) {
                    cval[0] = 1.0;
                    cycleval += g->xval[b_arc];
                } else {
                    cval[0] = -1.0;
                    rhs -= 1.0;
                    cycleval -= 1 - g->xval[b_arc];
                }
                                
                cind[1] = g->warc[in_edge[neighbor]];
                if (g->tail[in_edge[neighbor]] < g->head[in_edge[neighbor]]) {
                    cval[1] = 1.0;
                    cycleval += g->xval[in_edge[neighbor]];
                } else {
                    cval[1] = -1.0;
                    rhs -= 1.0;
                    cycleval -= 1 - g->xval[in_edge[neighbor]];
                }
                
                cind[2] = g->warc[arc];
                if (g->tail[arc] < g->head[arc]) {
                    cval[2] = 1.0;
                    cycleval += g->xval[arc];
                } else {
                    cval[2] = -1.0;
                    rhs -= 1.0;
                    cycleval -= 1 - g->xval[arc];
                }
                
                assert(cycleval < rhs + 1 + eps);

                if (cycleval > rhs + eps) {
                    for (i = 0; i < 3; i++) {
                        constraints[nconstr].cind[i] = cind[i];
                        constraints[nconstr].cval[i] = cval[i];
                    }
                    constraints[nconstr].rhs = rhs;
                    constraints[nconstr].violation = cycleval - rhs;
                    nconstr++;
                    if (nconstr == MAX_STORED_CONSTRAINTS) {
                        goto enough;
                    }
                }
            }

            arc = g->nextin[arc];
        }        
    }

enough:
    // shuffle constraints
    for (j=nconstr-1; j>0; j--) {
        i = gb_unif_rand(j);
        t = constraints[j];
        constraints[j] = constraints[i];
        constraints[i] = t;
    }

    int ngen = 0;
    int n = nconstr; if (n > max_gen_cycles) n = max_gen_cycles;
    for (i = 0; i < n; i++) {
        // insert row
        OsiRowCut rc;
        rc.setLb(-COIN_DBL_MAX);
        rc.setUb(constraints[i].rhs);
        rc.setRow(3, constraints[i].cind, constraints[i].cval);
        rc.setGloballyValid(true);
        osiCuts.insert(rc);
        ngen++;
    }
    
    if (pril) {
        if (ngen == nconstr)
            printf("%d 3-cycles\n", ngen);
        else
            printf("%d 3-cycles (of %d)\n", ngen, nconstr);
    }

    free(constraints);

    free(in_edge);

    free(b_edges);
    free(is_b_edge);

    free(current_neighborarc_stack);
    free(node_stack);
    free(nodepermutation);
    free(visited);

    num3cyclecuts += ngen;

    return ngen;
}

//
// HEURISTIC SEPARATION OF INTEGRAL CYCLES (LENGTH >=4)
// Method to shorten and add cycle contraint
//
//
int integralCycleSep_addBCycle(int b_arc, Graph* g, int* dfsPredArc, int minCycleLen, int maxCycleLen, OsiCuts& osiCuts) {
    
    int i, arc;
    int neighbor;

    int* partOfCycle = (int*) malloc(sizeof(int) * numbottomnodes); // whether node i is part of cycle found in DFS
    for (i = 0; i < numbottomnodes; i++)
        partOfCycle[i] = 0;

    int b_tail = g->tail[b_arc];
    int b_head = g->head[b_arc];

    int cycleLen = 1;
    partOfCycle[b_head] = 1;
    int node = b_tail;
    while (node != b_head) {
        partOfCycle[node] = 1;

        assert(dfsPredArc[node] >= 0);
        node = g->tail[dfsPredArc[node]];
        cycleLen++;
    }

    assert(cycleLen >= 3);
    #ifndef NDEBUG
    {
        int count = 0;
        for (int i = 0; i < numbottomnodes; i++)
            count += partOfCycle[i];
        assert(count == cycleLen);
    }
    #endif

    // find shortest path from b_head to b_tail using BFS, only using nodes that are part of original cycle
    int* bfsPredArc = (int*)malloc(numbottomnodes * sizeof(int)); // bfs predecessor
    int* bfsVisited = (int*)malloc(numbottomnodes * sizeof(int)); // to keep track of visited nodes
    for (i = 0; i < numbottomnodes; i++) {
        bfsPredArc[i] = -1;
        bfsVisited[i] = 0;
    }
    int* queue = (int*)malloc(cycleLen * sizeof(int)); // bfs queue
    int queueTail = 0;
    int queueHead = 0;

    // add b_head to queue
    queue[queueHead++] = b_head;
    bfsVisited[b_head] = 1;
    while (queueHead - queueTail > 0) { // while queue is not empty
        // pop tail of queue
        node = queue[queueTail++];
        assert(partOfCycle[node]);

        arc = g->firstout[node];
        // printf("node: %d\n", node);
        while (arc != -1) {
            neighbor = g->head[arc];

            if (partOfCycle[neighbor]) {
                // add neighbor to queue, if not visited
                if (bfsVisited[neighbor] == 0) {
                    bfsVisited[neighbor] = 1;
                    bfsPredArc[neighbor] = arc;
                    if (neighbor == b_tail) {
                        // found b_u
                        break;
                    }
                    assert(queueHead < cycleLen);
                    queue[queueHead++] = neighbor; // add to queue
                }
            }
            arc = g->nextout[arc];
        }
    }

    if (bfsPredArc[b_tail] == -1) {
        printf("ERROR, b_u not found");
        exit(1);
    }
    
    // reconstruct shortest path
    double rhs = 0.0;
    int newCycleLen = 0;
    int ncind = 0;
    node = b_tail;
    while (node != b_head) {
        if (!g->fixed[bfsPredArc[node]]) {
            cind[ncind] = g->warc[bfsPredArc[node]];
            if (g->tail[bfsPredArc[node]] < g->head[bfsPredArc[node]]) {
                cval[ncind] = 1.0;
            } else {
                cval[ncind] = -1.0;
                rhs -= 1.0;
            }
            ncind++;
        }
        newCycleLen++;

        node = g->tail[bfsPredArc[node]];
    }

    // add b-edge
    if (!g->fixed[b_arc]) {
        cind[ncind] = g->warc[b_arc];
        if (g->tail[b_arc] < g->head[b_arc]) {
            cval[ncind] = 1.0;
        } else {
            cval[ncind] = -1.0;
            rhs -= 1.0;
        }
        ncind++;
    }
    newCycleLen++;

    rhs += ncind - 1;

    // printf("Cycle new: %d\n", newCycleLen);

    assert(newCycleLen >= 3);
    if (newCycleLen < 3) {
        printf("Error: new cyclelen<3");
        exit(1);
    }

    // insert row
    int ngen = 0;

    if (newCycleLen >= minCycleLen && newCycleLen <= maxCycleLen) {
        OsiRowCut rc;
        rc.setLb(-COIN_DBL_MAX);
        rc.setUb(rhs);
        rc.setRow(ncind, cind, cval);
        rc.setGloballyValid(true);
        osiCuts.insert(rc);

        ngen++;
    }

    free(queue);
    free(bfsVisited);
    free(bfsPredArc);
    free(partOfCycle);

    return ngen;
}

//
// HEURISTIC SEPARATION OF INTEGRAL CYCLES (LENGTH >=4)
//
// todo, maybe only once in cutwidth track or stop early if no b cycle found
int integralCycleSep(Graph* g, OsiCuts& osiCuts, int minCycleSize, int numberIterations, int maxGen) {
    if (pril) {
        printf("INT-CYCLES ");fflush(stdout);
    }
    // DFS to find B-edges

    int i;
    int cycleLen, minCycleLen, maxCycleLen;
    int node, startNode;
    int current_neighborarc, neighbor;

    // arrays for dfs
    int* visited = (int*)malloc(sizeof(int) * numbottomnodes);
    int* predArc = (int*)malloc(sizeof(int) * numbottomnodes);

    int* nodepermutation = (int*)malloc(sizeof(int) * numbottomnodes);
    for (int i = 0; i < numbottomnodes; i++) {
        nodepermutation[i] = i;
    }

    // stacks for dfs
    int* node_stack = (int*) malloc(sizeof(int) * numbottomnodes);
    int* current_neighborarc_stack = (int*) malloc(sizeof(int) * numbottomnodes);
    int stack_size = 0;

    int ngen = 0;

    // shuffle before! just to be save
    shuffle(nodepermutation, numbottomnodes); // random node permutation

    for (cycleLen = minCycleSize; cycleLen < minCycleSize + numberIterations; cycleLen++) { // test smaller cycles first
        minCycleLen = cycleLen;
        maxCycleLen = cycleLen;
        if (cycleLen == minCycleSize + numberIterations - 1) // last iteration
            maxCycleLen = 10000000; // allow cycles of all lengths

        for (i = 0; i < numbottomnodes; i++) {
            visited[i] = 0;
            predArc[i] = -1;
        }

        for (int i = 0; i < numbottomnodes; i++) { // restart dfs on each node
            startNode = nodepermutation[i];
            if (visited[startNode] != 0)
                continue;
            node_stack[stack_size] = startNode; // add start node to stack
            current_neighborarc_stack[stack_size++] = g->firstout[startNode];

            while (stack_size > 0) {
                node = node_stack[stack_size - 1];
                if (visited[node] == 0) { // newly discovered
                    visited[node] = 1;
                }
                // traverse all neighbors
                current_neighborarc = current_neighborarc_stack[stack_size - 1];
                if (current_neighborarc != -1) {
                    current_neighborarc_stack[stack_size - 1] = g->nextout[current_neighborarc];
                    neighbor = g->head[current_neighborarc];
                    if (visited[neighbor] == 1) {
                        // cycle found
                        ngen += integralCycleSep_addBCycle(current_neighborarc, g, predArc, minCycleLen, maxCycleLen, osiCuts);
                    } else if (visited[neighbor] == 0) {
                        // add to stack
                        predArc[neighbor] = current_neighborarc;
                        node_stack[stack_size] = neighbor;
                        current_neighborarc_stack[stack_size++] = g->firstout[neighbor];
                    }
                } else {
                    visited[node] = 2;
                    stack_size--;
                }
                if (stack_size > numbottomnodes) {
                    printf("stack size %d>=nnod (%d)\n", stack_size, numbottomnodes);
                    exit(1);
                }

                if (ngen > maxGen) {
                    goto enough;
                }
            }
        }
        // if (ngen == 0) {
        //     if (pril) printf("No %d cycles\n", minCycleLen);
        // }
    }

enough:

    if (pril)
        printf("%d cycles\n", ngen);

    free(current_neighborarc_stack);
    free(node_stack);

    free(nodepermutation);
    free(predArc);
    free(visited);

    numlargecyclecuts += ngen;

    return ngen;
}

void buildSeparationGraph(Graph* g, const double* wxval, int onlyintegral) {
    int a, head, tail;

    // reset graph
    init_graph(*g, numbottomnodes);

    // insert edges
    for (a = 0; a < wnvar; a++) {
        head = whead[a];
        tail = wtail[a];
        double xval = fmin(fmax(wxval[a], 0), 1);
        // printf("%d -> %d\n", tail, head);
        if (xval > omeps || (!onlyintegral && xval > eps)) {
            // insert edge tail->head, xval
            insert_wedge(*g, tail, head, a, xval);
        }
        if (xval < eps || (!onlyintegral && xval < omeps)) {
            // insert edge head->tail, 1-xval
            insert_wedge(*g, head, tail, a, 1 - xval);
        }
    }
}

//
// CYCLE SEPARATION
//
void separate(const double* wxval, OsiCuts& osiCuts, double localLowerboundCbc, const CglTreeInfo info) {
    int a, i;
    double localLowerbound = localLowerboundCbc + objconst;
    currentCompLocalLowerbound = localLowerbound;
    numlps++;
    
    // tailing off
    const int MAX_ITER_SMALL_CHANGE = 18;
    if (!info.inTree) {
        if (fabs(localLowerbound - lastRootLowerbound) < 0.04) {
            numRootSmallObjChange++;
        } else {
            numRootSmallObjChange = 0;
        }
        if (numRootSmallObjChange > MAX_ITER_SMALL_CHANGE) {
            // only continue if integer feasible
            for (int i = 0; i < wnvar; i++) {
                if (wxval[i] > eps && wxval[i] < omeps) {
                    // branch
                    if (pril) printf("\nBranching at lowerbound %.4f\n\n", localLowerbound);
                    return;
                }
            }
        }
    }
    
    if (pril) {
        if (!info.inTree && numRootSmallObjChange)
            printf("[%.4lf,%.4lf] (small improvement since %d iterations) (time: %.0fs)\n",localLowerbound, activityrank_num_crossings, numRootSmallObjChange, processtime());
        else
            printf("[%.4lf,%.4lf] (time: %.0fs)\n",localLowerbound, activityrank_num_crossings, processtime());
    }
    
    // don't stay in same b&b node for long time
    if (info.inTree && info.pass > 5) {
        // only continue if integer feasible
        for (i = 0; i < wnvar; i++) {
            if (wxval[i] > eps && wxval[i] < omeps) {
                // branch
                return;
            }
        }
    }

    int ncuts = 0;

    // build separation graph including fractional edges
    buildSeparationGraph(sepgraph, wxval, false);

    // search for violated >=3 cycles in sepgraph
    ncuts += fast3CycleSep(sepgraph, osiCuts);

    if (ncuts < max_gen_cycles / 10) {
        // build graph for integral separation based on wgraph
        buildSeparationGraph(sepgraphint, wxval, true);

        // search for integral >=4 cycles
        ncuts += integralCycleSep(sepgraphint, osiCuts, 4, 3, max_gen_cycles); // todo: how many rounds
    }

    if (ncuts < max_gen_cycles / 4) {       
        // search for fractional cycles in sepgraph 
        ncuts += fascycle(sepgraph, true, osiCuts);
    }

    // no cuts found, search in combined graph
    if (!ncuts) {
        // insert fixed edges into sepgraphint
        for (a = 0; a < numfixed; a++) {
            insert_fixededge(*sepgraphint, ftail[a], fhead[a]);
        }
        // heuristic separation on combined graph
        // search for integral >=3 cycles
        ncuts += integralCycleSep(sepgraphint, osiCuts, 3, 1, 10);
        if (ncuts && pril) {
            printf("dfs combined has found cuts\n");
        }

        if (!ncuts) {
            int numberFractional = 0;
            for (i = 0; i < sepgraph->nedges; i++) {
                if (sepgraph->xval[i] > eps && sepgraph->xval[i] < omeps)
                    numberFractional++;
            }
            if (pril && numberFractional > 0) {
                printf("FRACTIONAL IN COMBINED: %d\n", numberFractional);
            }

            if (numberFractional > 0) {
                // insert fixed edges into sepgraph
                for (a = 0; a < numfixed; a++) {
                    insert_fixededge(*sepgraph, ftail[a], fhead[a]);
                }
                ncuts += fascycle(sepgraph, true, osiCuts); // run only on fractional edges of cgraph
                if (ncuts && pril) {
                    printf("fas combined has found cuts\n");
                }
            }
        }
    }

    if (ncuts == 0 || numRootSmallObjChange >= MAX_ITER_SMALL_CHANGE - 5) {
        separatemoebius = true;
    }

    // don't call on cutwidth instances to save time
    if (!iscutwidth && separatemoebius && (ncuts <= max_gen_cycles / 3 || numRootSmallObjChange == MAX_ITER_SMALL_CHANGE - 3)) {
        // moebius sep
        ncuts += moebiusSep(wxval, osiCuts);
    }

    lastRootLowerbound = localLowerbound;
}

//
// Branch & Cut
//

int bac()
{
    int a;   // arc variable
        
    double *lb;        // lower bounds
    double *ub;        // upper bounds
    double *startsol;  // starting solution vector for the b&c phase

    compstarttime = processtime();
    
    sepgraph = (Graph*) malloc(sizeof(Graph));
    allocate_graph(*sepgraph, numbottomnodes, wnvar * 2, numfixed);
    sepgraphint = (Graph*) malloc(sizeof(Graph));
    allocate_graph(*sepgraphint, numbottomnodes, wnvar, numfixed);

    lb       = (double *) malloc( wnvar * sizeof(double));
    ub       = (double *) malloc( wnvar * sizeof(double));
    startsol = (double *) malloc( wnvar * sizeof(double));

    cind      = (int *)    malloc( wnvar * sizeof(int));
    graph_cind= (int *)    malloc( wnvar * sizeof(int));
    cval      = (double *) malloc( wnvar * sizeof(double));
    
    separatemoebius = false;
    numRootSmallObjChange = 0;
    lastRootLowerbound  = 10e9;

    if (pril>1) {
        printf("\nSetting up cbc...\n\n");
        printf("Number of variables: %d\n",wnvar);
    }
    
    for (a=0; a<wnvar; a++) {
        lb[a] = 0.0;
        ub[a] = 1.0;
        if (rank[wtail[a]]<rank[whead[a]]) startsol[a] = 1; else startsol[a] = 0;
    }

    int generatorStats = true;
#ifdef SUBMISSION
    generatorStats = false;
#endif

    // initialize cbc model
    heuristic_callback heuristic = nullptr; // only use heurisic on exact track
    if (!iscutwidth) {
        heuristic = asexploitlp;
    }
    CbcInterface model(wnvar, lb, ub, obj, startsol, separate, heuristic, pril, generatorStats);
    CbcModel* model_ = model.model_;

    if (pril>1) printf("Model initialized.\n");
    
    cbcinittime += (processtime() - compstarttime);

#ifndef SUBMISSION
    remainingtime = programstarttime + processtimelimit - processtime(); // remaining time for MIP
    model_->setDblParam(CbcModel::CbcMaximumSeconds, remainingtime);
#endif

    model_->branchAndBound(pril > 1 ? 3 : 0);

    int optimumCrosscount = (int)(model_->getObjValue() + objconst + eps);
    if (pril>1) printf("\nOptimum solution: %d \n\n",optimumCrosscount);
    
    numbbnodes += model_->getNodeCount();

    if (generatorStats) {
        if (pril>1) printf("Statistics:\n");
        for (int iGenerator = 0; iGenerator < model_->numberCutGenerators(); iGenerator++) {
            CbcCutGenerator *generator = model_->cutGenerator(iGenerator);
            if (pril>1) printf("%s was tried %d times and created %d cuts of which %d were active after adding rounds of cuts",
                generator->cutGeneratorName(), generator->numberTimesEntered(), generator->numberCutsInTotal(),generator->numberCutsActive());
            if (generator->timing()) {
                if (pril>1) printf(" (%.2f seconds)\n", generator->timeInCutGenerator());
                if (std::string(generator->cutGeneratorName()).compare("LazyCuts") == 0) {
                    separationtimecallback += generator->timeInCutGenerator();
                } else {
                    separationtimeother += generator->timeInCutGenerator();
                }
            }
            else if (pril > 1)
                printf("\n");
        }
    }

    int returnCode = 0;

    if (model_->isProvenOptimal()) {
        // retrieve optimum solution
        const double* wxval = model_->getColSolution();

        // translate to rank
        double d;
        asexploitlp(wxval, false, d, nullptr);

        // check whether the exploiter did its thing     
        if ((int) (activityrank_num_crossings+eps) != optimumCrosscount) {
            if (pril) printf("Error: Cbc solution value does not match exploited value, trying again...\n");
            returnCode = 3; // bug in cbc, try to start again
        }

        //for (i=nnod-1; i>=0; i--) printf(" %d\n",activity[i]+numtopnodes+1);
    } else {
        timeout = true;
        returnCode = 1;
        #ifdef SUBMISSION
            // if solution is not proven optimal, quit immediately
            exit(1);
        #endif
    }
    
    bactime += (processtime() - compstarttime);
    
    free(lb);
    free(ub);
    free(cind);
    free(graph_cind);
    free(cval);
    free(startsol);

    free_graph(*sepgraphint);
    free(sepgraphint);
    free_graph(*sepgraph);
    free(sepgraph);

    return returnCode;
}

//
// MAIN PRGRAM
//

int main(int argc, char *argv[])
{
    // variables
    
    char inputline[1000]; // a line of the input
    
    int helpwanted     = false;  // true iff the -h command line option is applied
    int printsummary   = false;  // true iff a summary line should be printed
    int onlysummary    = false;  // true iff nothing but a summary line should be printed
    int fullcrosstable = false;  // true iff a full crosstable is constructed (not necessary, just a safety measure)
    int rseed          = 4711;   // random seed

    int comment;                // indicates a comment line
    int cnt;                    // return value of scanf
    int foo;                    // dummy
    int i,j,k;                  // running indices
    int v,w;                    // nodes on bottom layer
    int e;                      // edge variable for input graph
    int wa;                     // arc variable for weighted digraph
    int fa;                     // arc variable for fixed digraph
    int ev,ew;                  // edges incident with bottom nodes v and w
    int nv,nw;                  // top node neighbors of bottom nodes v and w
    int tne;                    // top node of edge e
    int bne;                    // bottom node of edge e
    int crvw,crwv;              // number of crossings between v and w if v before w and v after w, respectively
    int oni,onj;                // original names of bottom nodes i and j
    int mtonj;                  // mintopnode[onj]
    int mintop;                 // minimum and maximum top nodes in a component
    int maxtop;                 // minimum and maximum top nodes in a component
    int mintopj,maxtopj;        // minimum and maximum top nodes of node j
    int minbotindex;            // minimum index of component
    int numzerodegreenodes;     // number of zero degree nodes
    int comp;                   // component running index
    int ra;                     // rank variable
    int optbottomsequenceindex; // index in optbottomsequence
    int optnumcrossings;        // optimum number of crossings
    int maxtopv,mintopw;        // max/min top nodes of v/w

    double diff;  // difference of crvw and crwv

    int *topnode;           // array of top nodes
    int *bottomnode;        // array of bottom nodes
    int *firstedge;         // firstedge[v] is the first edge incident to bottom node v
    int *nextedge;          // nextedge[e] implements a linked list of edges
    int *mintopnode;        // minimum top node of a bottom node
    int *maxtopnode;        // maximum top node of a bottom node
    int *origname;          // original name of a bottom node
    int *component;         // component a bottom node belongs to
    int *compcard;          // component cardinalities
    int *degree;            // bottom node degrees
    int *optbottomsequence; // crossing minimum sequence for bottom nodes
    
    programstarttime = processtime();
    
    pril = 0;          // output (print) level
            
    #ifndef SUBMISSION
    // change values according to the command line options
    while (--argc) {
        if      (sscanf(argv[argc],"-p%d",&pril)==1);
        else if (sscanf(argv[argc],"-maxgen%d",&max_gen_cycles)==1);
        else if (sscanf(argv[argc],"-timelimit%lf",&processtimelimit)==1);
        else if (sscanf(argv[argc],"-rseed%d",&rseed)==1);
        else if (strcmp(argv[argc],"-onlysummary")==0) onlysummary = true;
        else if (strcmp(argv[argc],"-nooutput")==0) nooutput = true;
        else if (strcmp(argv[argc],"-fullcrosstable")==0) fullcrosstable = true;
        else if (strcmp(argv[argc],"-h")==0) helpwanted = true;
        else {
            printf("Usage: %s [-p<int>] [-maxgen<int>] [-timelimit<double>] [-rseed<int>] [-onlysummary] [-nooutput] [-fullcrosstable] [-h]\n",argv[0]);
            return 0;
        }
    }
    #endif

    gb_init_rand(rseed);
    
    if (helpwanted) {
        printf("\n\nUsage: %s [-p<int>] [-maxgen<int>] [-timelimit<int>] [-onlysummary] [-nooutput] [-fullcrosstable] [-h]\n\n",argv[0]);
        printf("   -p<int>             print level (default 0)\n");
        printf("   -maxgen<int>        maximum number of 3-cycles allowed in a separation step (default 500)\n");
        printf("   -timelimit<double>  time limit for entire computation in seconds (default no limit)\n");
        printf("   -rseed<int>         random seed (default 4711)\n");
        printf("   -onlysummary        print nothing but a summary line\n");
        printf("   -nooutput           don't write solution to stdout (default do unless p>0 or onlysummary)\n");
        printf("   -fullcrosstable     use full crosstable for verification (default don't)\n");
        printf("   -h                  get this help\n\n");
        return 0;
    }

    if ((pril>0)||(onlysummary)) {
        nooutput = true;
        printsummary = true;
    }
    
    starttime = processtime();
    
    // input of numtopnodes numbottomnodes numedges
    
    comment = true;
    do {
        inputline[0] = 0;
        cnt = scanf("%999[^\n]",inputline);
        if (cnt == EOF) {
            printf("unexpected eof\n");
            exit(2);
        }
        // Consume \n if next stdin char is a \n
        else cnt = scanf("%*1[\n]");
        // Use inputline;
        if (inputline[0]=='c') { // a comment line
            if (pril>2) printf("%s\n",inputline);
        }
        else if (inputline[0]=='p') { // the "p" line
            if (pril>1) printf("%s\n",inputline);
            comment = false;
            cnt = sscanf(inputline+5,"%d %d %d %d",&numtopnodes,&numbottomnodes,&numedges,&foo);
            if (pril > 1) {
                printf("\n%d numbers read in p-line.\n",cnt);
                printf("numtopnodes = %d\n",numtopnodes);
                printf("numbottomnodes = %d\n",numbottomnodes);
                printf("numedges = %d\n",numedges);
            }
            if (cnt==4) {
                // ignore the next numtopnodes+numbottomnodes numbers
                for (i=0; i<numtopnodes+numbottomnodes; i++)
                    cnt = scanf("%d",&foo);
                iscutwidth = true;
            } else {
                iscutwidth = false;
            }
        }
        else {
            printf("!!! corrupt input file !!!\n");
            exit(3);
        }
    } while (comment);

    // always apply fixing on cutwidth instances -> don't use initial heuristic
    if (!iscutwidth) {
        if (numbottomnodes <= 400)
            dofixing = false;
    }

    if (pril && !dofixing)
        printf("Fixing switched off\n");
    
    // allocations
    
    topnode           = (int *) malloc(numedges * sizeof(int));
    bottomnode        = (int *) malloc(numedges * sizeof(int));
    nextedge          = (int *) malloc(numedges * sizeof(int));
    firstedge         = (int *) malloc(numbottomnodes * sizeof(int));
    mintopnode        = (int *) malloc(numbottomnodes * sizeof(int));
    maxtopnode        = (int *) malloc(numbottomnodes * sizeof(int));
    origname          = (int *) malloc(numbottomnodes * sizeof(int));
    component         = (int *) malloc(numbottomnodes * sizeof(int));
    degree            = (int *) malloc(numbottomnodes * sizeof(int));
    compcard          = (int *) malloc(numbottomnodes * sizeof(int));
    node              = (int *) malloc(numbottomnodes * sizeof(int));
    poscct            = (int *) malloc(numbottomnodes * sizeof(int));
    optbottomsequence = (int *) malloc(numbottomnodes * sizeof(int));
    rank              = (int *) malloc(numbottomnodes*sizeof(int));
    activity          = (int *) malloc(numbottomnodes*sizeof(int));
    wfirst            = (int *) malloc(numbottomnodes*sizeof(int));
    wfirstin          = (int *) malloc(numbottomnodes*sizeof(int));
    ffirst            = (int *) malloc(numbottomnodes*sizeof(int));

    // construct 2-layer graph
    
    for (v=0; v<numbottomnodes; v++) {
        firstedge[v] = -1;
        mintopnode[v] = numtopnodes;
        maxtopnode[v] = 0;
        degree[v] = 0;
    }
    for (e=0; e<numedges; e++) {
        cnt = scanf("%d %d",&tne,&bne);
        // change xindexing to 0,1,2,... for both layers
        if (tne<bne) {
            topnode[e] = tne - 1;
            bottomnode[e] = bne - numtopnodes - 1;
        }
        else {
            topnode[e] = bne - 1;
            bottomnode[e] = tne - numtopnodes - 1;
        }
        tne = topnode[e];
        bne = bottomnode[e];
        nextedge[e] = firstedge[bne];
        firstedge[bne] = e;
        degree[bne]++;
        if (mintopnode[bne]>tne) mintopnode[bne] = tne;
        if (maxtopnode[bne]<tne) maxtopnode[bne] = tne;
    }
    
    // count zero degree nodes and assign large mintopnode and -1 component
    
    numzerodegreenodes = 0;
    for (i=0; i<numbottomnodes; i++) {
        if (!degree[i]) {
            numzerodegreenodes++;
            mintopnode[i] = numtopnodes;
            component[i] = -1;
        }
    }

    // show zero degree nodes
    
    if ((pril>1)&&(numzerodegreenodes)) {
        printf("\n\nThere are %d zero degree nodes:",numzerodegreenodes);
        for (i=0; i<numbottomnodes; i++) if (!degree[i]) printf(" %d",i);
        printf("\n");
    }
        
    // sort the bottom nodes according to mintopnode
    
    for (i=0; i<numbottomnodes; i++) origname[i] = i;

    for (j=1; j<numbottomnodes; j++) {
        onj = origname[j];
        mtonj = mintopnode[onj];
        i = j - 1;
        while ((i>=0)&&(mintopnode[origname[i]]>mtonj)) {
            origname[i+1] = origname[i];
            i--;
        }
        origname[i+1] = onj;
    }

    // find components
    
    numcomponents = 0;
    i = 0;
    oni = origname[i];
    minbotindex = 0;
    mintop = mintopnode[oni];
    maxtop = maxtopnode[oni];
    for (j=i+1; j<numbottomnodes-numzerodegreenodes; j++) {
        onj = origname[j];
        maxtopj = maxtopnode[onj];
        mintopj = mintopnode[onj];
        if (mintopj<maxtop) {
            // if (pril>2) printf("merging %d and %d\n",oni,onj);
            if (maxtopj>maxtop) maxtop = maxtopj;
        }
        else {
            // if (pril>2) {
            //     printf("\ncomponent %d with top interval [%d,%d] and bottom nodes { ",numcomponents,mintop,maxtop);
            //     for (k=minbotindex; k<j; k++) printf("%d ",origname[k]);
            //     printf("}\n");
            // }
            compcard[numcomponents] = j - minbotindex;
            for (k=minbotindex; k<j; k++) component[origname[k]] = numcomponents;
            minbotindex = j;
            mintop = mintopj;
            maxtop = maxtopj;
            numcomponents++;
        }
        oni = onj;
    }
    
    if (pril>2) {
        printf("\ncomponent %d with top interval [%d,%d] and bottom nodes { ",numcomponents,mintop,maxtop);
        for (k=minbotindex; k<numbottomnodes-numzerodegreenodes; k++) printf("%d ",origname[k]);
        printf("}\n");
    }

    compcard[numcomponents] = numbottomnodes - numzerodegreenodes - minbotindex;
    for (k=minbotindex; k<numbottomnodes-numzerodegreenodes; k++) component[origname[k]] = numcomponents;
    numcomponents++;
    
    if (pril>2) for (i=0; i<numbottomnodes; i++) printf("component[%d] = %d\n",origname[i],component[origname[i]]);
    
    if (pril>2) for (i=0; i<numcomponents; i++) printf("component %d has %d nodes\n",i,compcard[i]);

    if (pril) {
        if (numcomponents>1) printf("\nThere are %d components.\n\n",numcomponents);
        else printf("\nThere is %d component.\n\n",numcomponents);
    }
    
    free(origname);

    // make crosstable
    
    if (fullcrosstable) {
        
        crosstable = (int *) malloc(numbottomnodes*numbottomnodes*sizeof(int));
        
        if (crosstable==NULL) {
            printf("crosstable too large\n");
            exit(1000);
        }
        
        for (k=0; k<numbottomnodes*numbottomnodes; k++) crosstable[k] = 0;
        for (v=0; v<numbottomnodes-1; v++) {
            ev = firstedge[v];
            while (ev>=0) { // for all neighbours nv of v on top layer
                nv = topnode[ev];
                ev = nextedge[ev];
                for (w=v+1; w<numbottomnodes; w++) {
                    ew = firstedge[w];
                    while (ew>=0) { // for all neighbours nw of w on top layer
                        nw = topnode[ew];
                        ew = nextedge[ew];
                        if (nv!=nw) {
                            if (nv>nw) {
                                crosstable[v*numbottomnodes+w] += 1;
                            }
                            else {
                                crosstable[w*numbottomnodes+v] += 1;
                            }
                        } // if (nv!=nw)
                    } // while (ew>=0)
                } // for w
            } // while (ev>=0)
        } // for v
        
        // if (pril>2) {
        //     // show crosstable
        //     printf("\nCrosstable\n\n");
        //     printf("    ");
        //     for (v=0; v<numbottomnodes; v++) printf("%4d",v);
        //     printf("\n\n");
        //     for (v=0; v<numbottomnodes; v++) {
        //         printf("%4d",v);
        //         for (w=0; w<numbottomnodes; w++) {
        //             printf("%4d",crosstable[v*numbottomnodes+w]);
        //         }
        //         printf("\n\n");
        //     }
        // }
    }
    
    inputtranstime = processtime() - starttime;
    
    if (pril>1) for (comp=0; comp<numcomponents; comp++) printf("Component %d has %d nodes\n",comp,compcard[comp]);

    // solve all components
    
    optbottomsequenceindex = 0;
    optnumcrossings = 0;
    
    for (comp=0; comp<numcomponents; comp++) {
        
        int cbc_bug = false; // whether a bug in cbc occurred

        starttime = processtime();

        // initialize rank&activity for this component
        
        ra = 0;
        for (v=0; v<numbottomnodes; v++) if (component[v]==comp) {
            node[ra] = v;
            poscct[v] = ra;
            activity[ra] = v;
            rank[v] = ra;
            ra++;
        }

        nnod = compcard[comp];
        
        // make crosstable for this component
        
        compcrosstable = (int *) malloc(nnod*nnod*sizeof(int));
        
        for (k=0; k<nnod*nnod; k++) compcrosstable[k] = 0;
        for (i=0; i<nnod-1; i++) {
            v = node[i];
            maxtopv = maxtopnode[v];
            for (j=i+1; j<nnod; j++) {
                    w = node[j];
                    mintopw = mintopnode[w];

                    if (maxtopv>mintopw || !dofixing) { // check whether edge should be fixed
                        ev = firstedge[v];
                        while (ev>=0) { // for all neighbours nv of v on top layer
                            nv = topnode[ev];
                            ew = firstedge[w];
                            while (ew>=0) { // for all neighbours nw of w on top layer
                                nw = topnode[ew];
                                if (nv!=nw) {
                                    if (nv>nw) compcrosstable[i*nnod+j]++;
                                    else compcrosstable[j*nnod+i]++;
                                } // if (nv!=nw)
                                ew = nextedge[ew];
                            } // while (ew>=0)
                            ev = nextedge[ev];
                        } // while (ev>=0)
                    } // if maxtopv>mintopw
                    else {
                        // fix variable
                        if (degree[v] > 1 || degree[w] > 1 || maxtopv != mintopw) // don't fix if 0 crossings
                            compcrosstable[j*nnod+i] = 1; // arbitrary positive value
                    }
                } // for j
        } // for i

        // if (pril>2) {
        //     // show crosstable for this component
        //     printf("\nCrosstable for this component\n\n");
        //     printf("    ");
        //     for (i=0; i<nnod; i++) printf("%4d",node[i]);
        //     printf("\n\n");
        //     for (i=0; i<nnod; i++) {
        //         printf("%4d",node[i]);
        //         for (j=0; j<nnod; j++) {
        //             printf("%4d",compcrosstable[i*nnod+j]);
        //         }
        //         printf("\n\n");
        //     }
        // }

        componentcrosstabletime += processtime() - starttime;
        starttime = processtime();

        // determine maximum number of variables
        wnvar = 0;
        numfixed = 0;
        for (i=0; i<nnod-1; i++) {
            for (j=i+1; j<nnod; j++) {
                crvw = compcrosstable[i*nnod+j];
                crwv = compcrosstable[j*nnod+i];
                if (crvw!=crwv || !doarbitrary) {
                    if ((crvw==0||crwv==0) && crvw != crwv && dofixing) numfixed++;
                    else wnvar++;
                }
            } /* for w */
        } /* for v */

        // allocate digraphs

        wtail  = (unsigned short *)    malloc(wnvar*sizeof(unsigned short));
        whead  = (unsigned short *)    malloc(wnvar*sizeof(unsigned short));
        wnext  = (int *)    malloc(wnvar*sizeof(int));
        wnextin= (int *)    malloc(wnvar*sizeof(int));
        obj    = (double *) malloc(wnvar*sizeof(double));
        
        ftail  = (unsigned short *)    malloc(numfixed*sizeof(unsigned short));
        fhead  = (unsigned short *)    malloc(numfixed*sizeof(unsigned short));
        fnext  = (int *)    malloc(numfixed*sizeof(int));

        remainingtime = programstarttime + processtimelimit - processtime(); // remaining time
        
        // make digraphs
        for (v=0; v<numbottomnodes; v++) {
            wfirst[v] = -1;
            wfirstin[v] = -1;
            ffirst[v] = -1;
        }
        numequal = 0;
        wa = 0;
        fa = 0;
        objconst = 0.0;
        activityrank_num_crossings = 0;
        for (i=0; i<nnod-1; i++) {
            v = node[i];
            for (j=i+1; j<nnod; j++) {
                w = node[j];
                crvw = compcrosstable[i*nnod+j];
                crwv = compcrosstable[j*nnod+i];
                activityrank_num_crossings += crvw;
                if (crvw==crwv && doarbitrary) {
                    objconst += crvw;
                    numequal++;
                }
                else { // crvw!=crwv
                    if (!crvw && crvw != crwv && dofixing) {
                        fnext[fa] = ffirst[v];
                        ffirst[v] = fa;
                        ftail[fa] = v;
                        fhead[fa] = w;
                        
                        fa++;

                    } else if (!crwv && crvw != crwv && dofixing) {
                        fnext[fa] = ffirst[w];
                        ffirst[w] = fa;
                        ftail[fa] = w;
                        fhead[fa] = v;

                        fa++;

                    } else {
                        diff = crvw - crwv;
                        
                        wtail[wa] = v;
                        whead[wa] = w;
                        wnext[wa] = wfirst[v];
                        wfirst[v] = wa;

                        wnextin[wa] = wfirstin[w];
                        wfirstin[w] = wa;

                        obj[wa] = diff;

                        objconst += crwv;

                        wa++;
                    }
                } // crvw!=crwv
            } // for w
        } // for v
        

        numfixed = fa;
        wnvar = wa;

        if (pril) {
            printf("\n====== Component %d with %d nodes and %d arcs ======\n\n",comp,nnod,wnvar+numfixed+numequal);
            printf("%d variables eliminated due to arbitrary ordering\n",numequal);
            printf("%d variables eliminated due to fixed ordering\n",numfixed);
            printf("%d of %d variables remaining.\n\n",wnvar,wnvar+numfixed+numequal);
        }
        
    #ifdef TRALLALA
        // show weighted graph
        if (pril>2) {
            int a;
            printf("weighted graph\n");
            for (a=0; a<wnvar; a++) printf("obj[%d=(%d,%d)] = %6.2lf\n",a,wtail[a]+numtopnodes+1,whead[a]+numtopnodes+1,obj[a]);
            for (i=0; i<nnod; i++) {
                v = node[i];
                printf("neighbors of %d:",v+numtopnodes+1);
                a = wfirst[v];
                while (a!=-1) {
                    printf(" %d",whead[a]+numtopnodes+1);
                    a = wnext[a];
                }
                printf("\n");
            }
            
            printf("weightsum = %9.2lf\n",weightsum);
            printf("objconst = %9.2lf\n",objconst);
        }
        
        // show fix graph
        if (pril>2) {
            int a;
            printf("fix graph\n");
            for (a=0; a<numfixed; a++) printf("%d=(%d,%d)\n",a,ftail[a]+numtopnodes+1,fhead[a]+numtopnodes+1);
            for (i=0; i<nnod; i++) {
                v = node[i];
                printf("neighbors of %d:",v+numtopnodes+1);
                a = ffirst[v];
                while (a!=-1) {
                    printf(" %d",fhead[a]+numtopnodes+1);
                    a = fnext[a];
                }
                printf("\n");
            }
        }
        
        if (pril>2) {
            printf("\n");
            printf("Number of nodes:             %9d\n",nnod);
            printf("Number of arcs:              %9d\n",wnvar);
            printf("Objective function constant: %9.3lf\n",objconst);
            printf("Initial solution value:      %9.3lf\n",upperbound);
        }
    #endif

        componentgraphtime += processtime() - starttime;

        if (nnod==1) {
            if (pril>2) printf("one node component\n");
            optbottomsequence[optbottomsequenceindex++] = node[0];
        }
        else if (nnod==2) {
            v = node[0];
            w = node[1];
           
            if (compcrosstable[1] < compcrosstable[nnod]) {
                optbottomsequence[optbottomsequenceindex++] = v;
                optbottomsequence[optbottomsequenceindex++] = w;
                optnumcrossings += compcrosstable[1];
            } else {
                optbottomsequence[optbottomsequenceindex++] = w;
                optbottomsequence[optbottomsequenceindex++] = v;
                optnumcrossings += compcrosstable[nnod];
            }
        }

        else {
            // solve with branch and cut
            
            if (!dofixing && numcomponents==1 && wnvar > 10 /* to prevent bug in heuristic */) {
                // find good initial solution
                lopheu();
            }

            if (pril>1) printf("\ncalling bac\n");

            int returnCode = bac();
           
            if (returnCode == 1) { // timeout
                timeout = true;
                goto SUMMARY;
            }

            if (returnCode == 3) { // bug in Cbc, try everything again
                cbc_bug = true;
            } else {

                if (pril) printf("%.0lf crossings\n",activityrank_num_crossings);
                optnumcrossings += (int) floor(activityrank_num_crossings+eps);
                for (i=0; i<nnod; i++) optbottomsequence[optbottomsequenceindex++] = activity[i];

            }

            // if (pril>2) {
            //     for (i=0; i<nnod; i++) printf(" %d",activity[i]+numtopnodes+1);
            //     printf("\n");
            // }
        }
        free(compcrosstable);
        
        free(wtail);
        free(whead);
        free(wnext);
        free(wnextin);
        free(obj);
        
        free(ftail);
        free(fhead);
        free(fnext);

        if (cbc_bug) {
            comp--; // try everything again due to bug in Cbc
        }

    } // for (comp=0; comp<numcomponents; comp++)
    
    if (optbottomsequenceindex<numbottomnodes) {
        // patch up with zero degree bottom nodes
        if (pril>2) printf("\nzero degree nodes\n");

        for (v=0; v<numbottomnodes; v++) if (component[v]==-1) {
            optbottomsequence[optbottomsequenceindex++] = v;
            if (pril>2) printf(" %d",v+numtopnodes+1);
        }
        if (pril>2) printf("\n");
    }
    if (optbottomsequenceindex!=numbottomnodes) printf("!!! optbottomsequenceindex!=numbottomnodes !!!\n");
    
    // verbose output
    
    totaltime = processtime() - programstarttime;
    
    if (pril) {
        printf("\n======= Summary =======\n\n");
        printf("Number of top nodes            %10d\n",numtopnodes);
        printf("Number of bottom nodes         %10d\n",numbottomnodes);
        printf("Number of nodes                %10d\n",numtopnodes+numbottomnodes);
        printf("Number of edges                %10d\n",numedges);
        printf("Number of components           %10d\n",numcomponents);
        printf("Number of 3-cycle cuts         %10d\n",num3cyclecuts);
        printf("Number of >3-cycle cuts        %10d\n",numlargecyclecuts);
        printf("Number of frac cycle cuts      %10d\n",numfraccyclecuts);
        printf("Number of moebius cuts         %10d\n",nummoebiuscuts);
        printf("Number of LPs                  %10d\n",numlps);
        printf("Number of branch&bound nodes   %10d\n",((int)numbbnodes));
        printf("Number of crossings            %10d\n\n",optnumcrossings);
        
        printf("Input and transformation time         "); printtime(inputtranstime);
        printf(" ( %6.2lf%% )\n\n",(100.0*inputtranstime)/totaltime);
        printf("Component crosstable time             "); printtime(componentcrosstabletime);
        printf(" ( %6.2lf%% )\n",(100.0*componentcrosstabletime)/totaltime);
        printf("Component graph time                  "); printtime(componentgraphtime);
        printf(" ( %6.2lf%% )\n\n",(100.0*componentgraphtime)/totaltime);
        printf("Total branch&cut time                 "); printtime(bactime);
        printf(" ( %6.2lf%% )\n\n",(100.0*bactime)/totaltime);
        printf("  including\n");
        printf("    exploitation time          "); printtime(exploitationtime);
        printf(" (%6.2lf%% )\n",(100.0*exploitationtime)/bactime);
        printf("    callback separation time     "); printtime(separationtimecallback);
        printf(" (%6.2lf%% )\n",(100.0*separationtimecallback)/bactime);
        printf("    including\n");
        printf("        moebius separation time     "); printtime(separationtimemoebius);
        printf(" (%6.2lf%% )\n",(100.0*separationtimemoebius)/separationtimecallback);
        printf("    other separation time      "); printtime(separationtimeother);
        printf(" (%6.2lf%% )\n",(100.0*separationtimeother)/bactime);
        printf("    cbc init time              "); printtime(cbcinittime);
        printf(" (%6.2lf%% )\n\n",(100.0*cbcinittime)/bactime);
        printf("Total computation time                "); printtime(totaltime);
        printf(" ( 100.00%% )\n\n");
    }
    
    if (fullcrosstable) {
        // verify with full crosstable
        k = 0;
        for (i=0; i<numbottomnodes-1; i++) {
            v = optbottomsequence[i];
            for (j=i+1; j<numbottomnodes; j++) {
                w = optbottomsequence[j];
                k += crosstable[v*numbottomnodes+w];
            }
        }
        if (k!=optnumcrossings) printf("\n !!! Wrong sequence (%d) !!!\n\n",k);
    }
    
    // output solution
    if (!nooutput) {
        for (i=0; i<numbottomnodes; i++) printf("%d\n",optbottomsequence[i]+numtopnodes+1);
    }

    // print summary line
SUMMARY:
    totaltime = processtime() - programstarttime;
    if (printsummary) {
        if (!timeout) {
            printf("& %6d & %6d & %6d & %6d & %6d & %6d & %6d ", numlps, numcomponents, numbbnodes, num3cyclecuts, numlargecyclecuts, numfraccyclecuts, nummoebiuscuts);
            if (!timeout)
                printf("& \\textbf{%10d} ", optnumcrossings);
            else
                printf("& [%d, %d] ", (int)(currentCompLocalLowerbound+eps),(int)(activityrank_num_crossings+eps));
            printf(" & %6.2lf",(100.0*separationtimecallback)/totaltime);
            printf(" & %6.2lf",(100.0*separationtimemoebius)/totaltime);
            printf(" & %6.2lf &",(100.0*separationtimeother)/totaltime);
            printtime(totaltime);
            printf(" \\\\\n");

        } else {
            printf("& \\cbl{%6d} & \\cbl{%6d} & \\cbl{%6d} & \\cbl{%6d} & \\cbl{%6d} & \\cbl{%6d} & \\cbl{%6d} ", numlps, numcomponents, numbbnodes, num3cyclecuts, numlargecyclecuts, numfraccyclecuts, nummoebiuscuts);
            if (!timeout)
                printf("& \\cbl{%10d} ", optnumcrossings);
            else
                printf("& \\cbl{[%d, %d]} ", (int)(currentCompLocalLowerbound+eps),(int)(activityrank_num_crossings+eps));
            printf(" & \\cbl{%6.2lf}",(100.0*separationtimecallback)/totaltime);
            printf(" & \\cbl{%6.2lf}",(100.0*separationtimemoebius)/totaltime);
            printf(" & \\cbl{%6.2lf} &",(100.0*separationtimeother)/totaltime);
            printtime(totaltime);
            printf(" \\\\\n");
        }
    }
    
    free(topnode);
    free(bottomnode);
    free(nextedge);
    free(firstedge);
    free(wfirst);
    free(wfirstin);
    free(ffirst);
    free(rank);
    free(activity);
    if (fullcrosstable) free(crosstable);
    free(mintopnode);
    free(maxtopnode);
    free(component);
    free(compcard);
    free(node);
    free(poscct);
    free(optbottomsequence);
    free(degree);

} // main
