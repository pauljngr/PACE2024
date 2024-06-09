#include <CbcModel.hpp>
#include <CglTreeInfo.hpp>
#include <ClpSimplex.hpp>
#include <OsiCuts.hpp>

#include "CbcHeuristic.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglPreProcess.hpp"
#include "CglProbing.hpp"
#include "CglRedSplit.hpp"

typedef void (*lazy_constr_callback)(const double* xval, OsiCuts& osiCuts, double localLowerbound, const CglTreeInfo info);
typedef int (*heuristic_callback)(const double* xval, int cbcCallback, double& solutionValue, double* betterSolution);

// Class to disallow strong branching solutions
#include "CbcFeasibilityBase.hpp"
class CbcFeasibilityNoStrong : public CbcFeasibilityBase {
   public:
    // Default Constructor
    CbcFeasibilityNoStrong(){};

    virtual ~CbcFeasibilityNoStrong(){};
    // Copy constructor
    CbcFeasibilityNoStrong(const CbcFeasibilityNoStrong& rhs){};

    // Assignment operator
    CbcFeasibilityNoStrong& operator=(const CbcFeasibilityNoStrong& rhs) {
        return *this;
    };

    /// Clone
    virtual CbcFeasibilityBase* clone() const {
        return new CbcFeasibilityNoStrong();
    };

    /**
       On input mode:
       0 - called after a solve but before any cuts
       -1 - called after strong branching
       Returns :
       0 - no opinion
       -1 pretend infeasible
       1 pretend integer solution
    */
    virtual int feasible(CbcModel* model, int mode) {
        return mode;
    };
};

class ProblemHeuristic : public CbcHeuristic {
   public:
    ProblemHeuristic(CbcModel* model, heuristic_callback callback);

    virtual CbcHeuristic* clone() const;
    ProblemHeuristic(ProblemHeuristic& other);

    virtual int solution(double& solutionValue, double* betterSolution);

    virtual void setModel(CbcModel* model);
    virtual void resetModel(CbcModel* model);

    virtual bool shouldHeurRun(int whereFrom);

   private:
    CbcModel* model_;
    heuristic_callback heuristicCallback_;
};

class LazyCutGenerator : public CglCutGenerator {
   public:
    virtual void generateCuts(const OsiSolverInterface& si, OsiCuts& cs, const CglTreeInfo info = CglTreeInfo());
    LazyCutGenerator(lazy_constr_callback callback, CbcModel* model);

    virtual CglCutGenerator* clone() const;
    LazyCutGenerator(const LazyCutGenerator&);

   private:
    lazy_constr_callback cutCallback_;
    CbcModel* model_;
};

class CbcInterface {
   public:
    CbcInterface(int nvar, const double* lb, const double* ub, const double* obj, const double* initialSolution, lazy_constr_callback lazyConstrCallback, heuristic_callback heuristicCallback, int pril, int generatorTiming);

    ~CbcInterface();

    CbcModel* model_;

   private:
    bool initialized_ = false;
    OsiClpSolverInterface* solver_;

    // CglProbing* probing_;
    CglGomory* gomory_;
    CglKnapsackCover* knacksackCover_;
    CglRedSplit* redSplit_;
    CglClique* clique_;
    CglMixedIntegerRounding2* mixedInt_;
    CglFlowCover* flowCover_;

    LazyCutGenerator* lazyCuts_;
    ProblemHeuristic* heuristic_;

    OsiBabSolver* osiBab_;
    CbcFeasibilityNoStrong* noStrong_;
};
