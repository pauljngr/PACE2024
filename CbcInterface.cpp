//
// Utility functions for using Cbc api
//
#include "CbcInterface.hpp"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/resource.h>


// For Branch and bound
#include "CbcModel.hpp"
#include "OsiSolverInterface.hpp"
#include "OsiClpSolverInterface.hpp"
#include "CbcCutGenerator.hpp"
#include "CbcHeuristicLocal.hpp"
#include "CbcStrategy.hpp"
#include "CoinBuild.hpp"
#include "CoinPragma.hpp"

// Cuts
#include "CbcHeuristic.hpp"
#include "CglClique.hpp"
#include "CglFlowCover.hpp"
#include "CglGomory.hpp"
#include "CglKnapsackCover.hpp"
#include "CglMixedIntegerRounding2.hpp"
#include "CglPreProcess.hpp"
//#include "CglProbing.hpp"
#include "CglRedSplit.hpp"
#include "CglStored.hpp"
#include "CoinTime.hpp"
#include "OsiAuxInfo.hpp"

ProblemHeuristic::ProblemHeuristic(ProblemHeuristic& rhs) {
    this->heuristicCallback_ = rhs.heuristicCallback_;
}
CbcHeuristic* ProblemHeuristic::clone() const {
    ProblemHeuristic* gen = new ProblemHeuristic(this->model_, this->heuristicCallback_);
    return static_cast<CbcHeuristic*>(gen);
}
ProblemHeuristic::ProblemHeuristic(CbcModel* model, heuristic_callback callback)
    : model_(model), heuristicCallback_(callback) {}

bool ProblemHeuristic::shouldHeurRun(int whereFrom) {
    // return whereFrom & 2; // only after root
    // return 0;
    return 1;
}

int ProblemHeuristic::solution(double& solutionValue, double* betterSolution) {
    OsiSolverInterface* solver = model_->solver();
    const double* xval = solver->getColSolution();
    // printf("HEUR\n");
    return heuristicCallback_(xval, true, solutionValue, betterSolution);
}
void ProblemHeuristic::setModel(CbcModel* model) {
    this->model_ = model;
}
void ProblemHeuristic::resetModel(CbcModel* model) {
    this->model_ = model;
}

void LazyCutGenerator::generateCuts(const OsiSolverInterface& si, OsiCuts& cs, const CglTreeInfo info) {
    const double* xval = si.getColSolution();
    cutCallback_(xval, cs, si.getObjValue(), info);
}
LazyCutGenerator::LazyCutGenerator(lazy_constr_callback callback, CbcModel* model) : cutCallback_(callback), model_(model) {}

CglCutGenerator* LazyCutGenerator::clone() const {
    LazyCutGenerator* gen = new LazyCutGenerator(this->cutCallback_, this->model_);
    return static_cast<CglCutGenerator*>(gen);
}
LazyCutGenerator::LazyCutGenerator(const LazyCutGenerator& rhs) {
    this->cutCallback_ = rhs.cutCallback_;
}

void doBranchAndBound(CbcModel* model) {
    model->branchAndBound(3);
}

CbcInterface::CbcInterface(int nvar, const double* lb, const double* ub, const double* obj, const double* initialSolution, lazy_constr_callback lazyConstrCallback, heuristic_callback heuristicCallback, int pril, int generatorTiming) {
    solver_ = new OsiClpSolverInterface();

    // Build Cbc model
    int* starts = new int[nvar + 1];
    for (int i = 0; i < nvar + 1; i++)
        starts[i] = 0;
    int r = 0;
    double v = 0;
    solver_->addCols(nvar, starts, &r, &v, lb, ub, obj);
    delete[] starts;

    for (int i = 0; i < nvar; i++) {
        solver_->setInteger(i);
    }

    if (pril)
        printf("%d vars, %d constr\n", solver_->getNumCols(), solver_->getNumRows());

    model_ = new CbcModel(*solver_);
    if (pril > 1)
        model_->messageHandler()->setLogLevel(3);
    else
        model_->messageHandler()->setLogLevel(0);

    model_->solver()->messageHandler()->setLogLevel(0);

    // todo, probing takes too long
//    probing_ = new CglProbing();
//    probing_->setUsingObjective(true);
//    probing_->setMaxPass(1);
//    probing_->setMaxPassRoot(5);
//    // Number of unsatisfied variables to look at
//    probing_->setMaxProbe(10);
//    probing_->setMaxProbeRoot(1000);
//    // How far to follow the consequences
//    probing_->setMaxLook(50);
//    probing_->setMaxLookRoot(500);
//    // Only look at rows with fewer than this number of elements
//    probing_->setMaxElements(200);
//    probing_->setRowCuts(3);

    gomory_ = new CglGomory();
    // try larger limit
    gomory_->setLimit(300);

    knacksackCover_ = new CglKnapsackCover();

    redSplit_ = new CglRedSplit();
    // try larger limit
    redSplit_->setLimit(200);

    clique_ = new CglClique();
    clique_->setStarCliqueReport(false);
    clique_->setRowCliqueReport(false);

    mixedInt_ = new CglMixedIntegerRounding2();
    flowCover_ = new CglFlowCover();

    // Add in cut generators
    // Experiment with -1 and -99 etc
//    model_->addCutGenerator(probing_, -1, "Probing");
    model_->addCutGenerator(gomory_, 1, "Gomory");
    model_->addCutGenerator(knacksackCover_, -1, "Knapsack"); // todo? https://github.com/coin-or/Cbc.old/issues/116
    // model.addCutGenerator(&generator4,-1,"RedSplit");
    // model.addCutGenerator(&generator5,-1,"Clique");
    model_->addCutGenerator(clique_, 1, "Clique");
    model_->addCutGenerator(flowCover_, -1, "FlowCover");
    model_->addCutGenerator(mixedInt_, -1, "MixedIntegerRounding");

    lazyCuts_ = new LazyCutGenerator(lazyConstrCallback, model_);
    // Add lazy cuts (making sure at all depths)
    model_->addCutGenerator(lazyCuts_, 1, "LazyCuts", true, true, false, -100, 1, -1); // also call on solution!!!

    // Say we want timings
    int numberGenerators = model_->numberCutGenerators();
    if (generatorTiming) {
        for (int iGenerator = 0; iGenerator < numberGenerators; iGenerator++) {
            CbcCutGenerator* generator = model_->cutGenerator(iGenerator);
            generator->setTiming(true);
        }
    }
    // "don't allow dual stuff"
    OsiClpSolverInterface* osiclp = dynamic_cast<OsiClpSolverInterface*>(model_->solver());
    osiclp->setSpecialOptions(osiclp->specialOptions() | 262144);

    // Uncommenting this should switch off all CBC messages
    // model.messagesPointer()->setDetailMessages(10,10000,NULL);

    // provide initial solution
    if (initialSolution) {
        // compute objective
        double objval = 0;
        for (int i = 0; i < nvar; i++) {
            objval += obj[i] * initialSolution[i];
        }
        model_->setBestSolution(initialSolution, nvar, objval, false);
    }

    // add heuristic
    if (heuristicCallback) {
        heuristic_ = new ProblemHeuristic(model_, heuristicCallback);
        model_->addHeuristic(heuristic_, "ASEXPLOIT", 0);
    } else {
        heuristic_ = nullptr;
    }

    // Do initial solve to continuous
    model_->initialSolve();
    /*  You need the next few lines -
        a) so that cut generator will always be called again if it generated cuts
        b) it is known that matrix is not enough to define problem so do cuts even
        if it looks integer feasible at continuous optimum.
        c) a solution found by strong branching will be ignored.
        d) don't recompute a solution once found
    */
    // Make sure cut generator called correctly (a)
    int iGenerator = numberGenerators - 1;
    model_->cutGenerator(iGenerator)->setMustCallAgain(true);
    // Say cuts needed at continuous (b)
    osiBab_ = new OsiBabSolver();
    osiBab_->setSolverType(4);
    // owing to bug must set after initialSolve
    model_->passInSolverCharacteristics(osiBab_);

    // Say no to all solutions by strong branching (c)
    noStrong_ = new CbcFeasibilityNoStrong();
    model_->setProblemFeasibility(noStrong_);

    // Say don't recompute solution d)
    // model.setSpecialOptions(4);

    // Switch off strong branching if wanted
    model_->setNumberStrong(0);
    //model_->setNumberStrong(20);
    // model.setNumberStrong(5);
    //model_->setNumberBeforeTrust(5);

    // todo
    // model.solver()->setIntParam(OsiMaxNumIterationHotStart, 100);

    if (pril)
        printf("Initialized model\n");

    initialized_ = true;
}

CbcInterface::~CbcInterface() {
    if (!initialized_)
        return;
    delete noStrong_;
    delete osiBab_;
    if (heuristic_)
        delete heuristic_;
    delete lazyCuts_;

    delete flowCover_;
    delete mixedInt_;
    delete clique_;
    delete redSplit_;
    delete knacksackCover_;
    delete gomory_;
    //delete probing_;

    delete solver_;
    delete model_;
}
