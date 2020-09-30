//
//  placement.cpp
//  Evolutionary placement: adding taxa to a phylogenetic tree
//  This file created by James Barbetti on 25/8/20, but:
//  1. addNewTaxaToTree was formerly in phylotree.cpp;
//  2. addTaxonML likewise (and was the work of BUI Quang Minh)
//
#include "phylotree.h"

namespace {
    enum CostFunction {
        MAXIMUM_PARSIMONY,           //maximum parsimony  (each taxon, each possible insertion place)
        SANKOFF_PARSIMONY,           //ditto, but using Sankoff parsimony
        MAXIMUM_LIKELIHOOD_MIDPOINT, //maximum likelihood (at midpoint of existing branch)
        MAXIMUM_LIKELIHOOD_ANYWHERE  //maximum likelihood (anyhwere in existing branch)
    };
    std::string getIncrementalParameter(const char letter, const char* defaultValue) {
        const std::string& inc = Params::getInstance().incremental_method;
        std::string answer = defaultValue;
        int braceLevel = 0;
        int i;
        for (i=0; i<inc.length(); ++i) {
            if (inc[i]==letter && braceLevel==0) {
                break;
            } else if (inc[i]=='{') {
                ++braceLevel;
            } else if (inc[i]=='}') {
                --braceLevel;
            }
        }
        if (i==inc.length()) {
            return answer;  //Didn't find it
        }
        ++i;
        defaultValue = "";
        int j;
        for (j=i; j<inc.length(); ++j) {
            if (inc[j]=='+' && braceLevel==0) {
                break;
            } else if (inc[j]=='-' && braceLevel==0) {
                break;
            } else if (inc[j]=='{') {
                ++braceLevel;
            } else if (inc[j]=='}') {
                --braceLevel;
            }
        }
        answer = inc.substr(i, j-i);
        if (!answer.empty() && answer[0]=='{'
            && answer[answer.length()-1]=='}' ) {
            answer = answer.substr(1, answer.length()-2);
        }
        return answer;
    }
    size_t getIncrementalParameter(const char letter, size_t defaultValue) {
        auto s = getIncrementalParameter(letter, "");
        if (s.empty()) {
            return defaultValue;
        }
        int i = convert_int_nothrow(s.c_str(), defaultValue);
        if (i<0) {
            return defaultValue;
        }
        return static_cast<size_t>(i);
    }
    size_t getNumberOfTaxaToRemove(size_t countOfTaxa) {
        if (countOfTaxa<4) {
            return 0;
        }
        string removalString = getIncrementalParameter('R', "");
        size_t len = removalString.length();
        if (len==0) {
            return 0;
        }
        size_t numberToRemove;
        if (removalString[len-1] == '%') {
            removalString = removalString.substr(0, len-1);
            double percent = convert_double_nothrow(removalString.c_str(), 0);
            if (percent<100.0/countOfTaxa) {
                return 0;
            } else if (100.0<=percent) {
                return 0; //Just ignore it. Todo: warn it's being ignored.
            }
            numberToRemove = (size_t) floor(percent * countOfTaxa / 100.0 + .5 );
        } else {
            numberToRemove = convert_int_nothrow(removalString.c_str(),0);
        }
        if (numberToRemove<1 || countOfTaxa <= numberToRemove+3) {
            return 0;
        }
        return numberToRemove;
    }
    CostFunction getCostFunction() {
        auto cf = getIncrementalParameter('C', "MP");
        if (cf=="ML") {
            return MAXIMUM_LIKELIHOOD_MIDPOINT;
        } else if (cf=="FML") {
            return MAXIMUM_LIKELIHOOD_ANYWHERE;
        } else if (cf=="SMP") {
            return SANKOFF_PARSIMONY;
        }
        return MAXIMUM_PARSIMONY;
    }
    enum LocalCleanup {
        NO_LOCAL_CLEANUP
    };
    LocalCleanup getLocalCleanupAlgorithm() {
        auto f = getIncrementalParameter('L', "");
        return NO_LOCAL_CLEANUP;
    }
    size_t getTaxaPerBatch(size_t totalTaxa) {
        size_t taxaPerBatch = getIncrementalParameter('B', 1);
        if (taxaPerBatch==0) {
            taxaPerBatch = totalTaxa;
        }
        return taxaPerBatch;
    }
    size_t getInsertsPerBatch(size_t totalTaxa, size_t taxaPerBatch) {
        string insertString = getIncrementalParameter('I', "");
        size_t len = insertString.length();
        if (len==0) {
            return 0;
        }
        size_t numberToInsert;
        if (insertString[len-1] == '%') {
            insertString = insertString.substr(0, len-1);
            double percent = convert_double_nothrow(insertString.c_str(), 0);
            if (percent<100.0/taxaPerBatch) {
                return 1;
            } else if (100.0<=percent) {
                return taxaPerBatch; //Just ignore it. Todo: warn it's being ignored.
            }
            numberToInsert = (size_t) floor(percent * taxaPerBatch / 100.0 + .5 );
        } else {
            numberToInsert = convert_int_nothrow(insertString.c_str(),0);
        }
        if (numberToInsert < 1 ) {
            numberToInsert = taxaPerBatch;
        }
        return numberToInsert;
    }
    enum BatchCleanup {
        NO_BATCH_CLEANUP
    };
    BatchCleanup getBatchCleanupAlgorithm() {
        auto f = getIncrementalParameter('A', "");
        return NO_BATCH_CLEANUP;
    }
    enum GlobalCleanup {
        NO_GLOBAL_CLEANUP
    };
    GlobalCleanup getGlobalCleanupAlgorithm() {
        auto f = getIncrementalParameter('T', "");
        return NO_GLOBAL_CLEANUP;
    }
};

/****************************************************************************
 Stepwise addition (greedy) by maximum likelihood
 ****************************************************************************/

double PhyloTree::recomputeParsimonyBranchLength(PhyloNode* fromNode, PhyloNode* toNode) {
    PhyloNeighbor* nei     = fromNode->findNeighbor(toNode);
    PhyloNeighbor* backnei = toNode->findNeighbor(fromNode);
    int       branch_subst = 0;
    computeParsimonyBranchFast(nei, fromNode, &branch_subst);
    double uncorrected_length = (branch_subst > 0)
                    ? ((double) branch_subst / getAlnNSite())
                    : (1.0 / getAlnNSite());
    double alpha    = (site_rate) ? site_rate->getGammaShape() : 1.0;
    nei->length     = correctBranchLengthF81(uncorrected_length, alpha);
    backnei->length = nei->length;
    return nei->length;
}

double PhyloTree::addTaxonML(PhyloNode* added_taxon,     PhyloNode *added_node,
                             PhyloNode* node,            PhyloNode* dad,
                             bool isAddedAtMidpoint,
                             PhyloNode* &target_node,    PhyloNode* &target_dad,
                             double& len_to_new_taxon,
                             double& len_to_target_node, double& len_to_target_dad) {

    Neighbor *dad_nei = dad->findNeighbor(node);

    //link the new interior node into the middle of the branch node-dad:
    //
    //   dad <---*---> added_node <---*---> node
    //                      ^
    //                      |
    //                      V
    //                 added_taxon
    //
    double len     = dad_nei->length;
    double halfLen = 0.5 * len;
    node->updateNeighbor(dad, added_node, halfLen);
    dad->updateNeighbor(node, added_node, halfLen);
    added_node->updateNeighbor(DUMMY_NODE_1, node, halfLen);
    added_node->updateNeighbor(DUMMY_NODE_2, dad,  halfLen);
    added_node->updateNeighbor(added_taxon, added_taxon, -1);
    added_taxon->updateNeighbor(added_node, added_node, -1);
    
    LOG_LINE(VB_DEBUG, "  Placement branch length " << len);

    FOR_EACH_PHYLO_NEIGHBOR(added_node, nullptr, it, nei) {
        nei->clearComputedFlags();
        nei->getNode()->findNeighbor(added_node)->clearComputedFlags();
    }
    
    //compute the likelihood
    PhyloNeighbor* nei;
    double best_score = 0;
    if (isAddedAtMidpoint) {
        len_to_new_taxon   = recomputeParsimonyBranchLength(added_taxon, added_node);
        LOG_LINE(VB_DEBUG, "  Parsimony taxon->interior length " << len_to_new_taxon );
        nei                = added_taxon->findNeighbor(added_node);
        best_score         = computeLikelihoodBranch(nei, added_taxon);
        LOG_LINE(VB_DEBUG, "  Traversal info size is " << traversal_info.size());
        LOG_LINE(VB_DEBUG, "  Likelihood before optimization " << best_score);
        optimizeOneBranch(added_taxon, added_node, false, 20);
        len_to_target_dad  = halfLen;
        len_to_target_node = halfLen;
        len_to_new_taxon   = nei->length;
        best_score         = computeLikelihoodFromBuffer();
        LOG_LINE(VB_DEBUG, "  Likelihood after optimization " << best_score
                           << " (len = " << len_to_new_taxon << ")");
    }
    else {
        len_to_new_taxon = recomputeParsimonyBranchLength(added_taxon, added_node);
        optimizeOneBranch(added_node,  dad, false, 20);
        nei                = added_node->findNeighbor(dad);
        len_to_target_dad  = nei->length;

        optimizeOneBranch(added_node,  node, false, 20);
        nei                = added_node->findNeighbor(node);
        len_to_target_node = nei->length;
        
        optimizeOneBranch(added_taxon,  added_node, false, 20);
        nei                = added_node->findNeighbor(added_taxon);
        best_score         = computeLikelihoodFromBuffer();
        len_to_new_taxon   = nei->length;
    }
    target_node        = node;
    target_dad         = dad;
    LOG_LINE(VB_DEBUG, "  ML Lengths " << len_to_target_dad
             << ", " << len_to_target_node << ", " << len_to_new_taxon << std::endl);

    //unlink the added node
    node->updateNeighbor(added_node, dad, len);
    dad->updateNeighbor(added_node, node, len);
    added_node->updateNeighbor(node, DUMMY_NODE_1, halfLen);
    added_node->updateNeighbor(dad,  DUMMY_NODE_2, halfLen);
    node->findNeighbor(dad)->clearComputedFlags();
    dad->findNeighbor(node)->clearComputedFlags();
    trackProgress(1.0);

    //now traverse the tree downwards
    FOR_EACH_ADJACENT_PHYLO_NODE(node, dad, it, child) {
        PhyloNode* target_node2 = nullptr;
        PhyloNode* target_dad2  = nullptr;
        double     len_child    = 0;
        double     len_node     = 0;
        double     len_dad      = 0;
        double score = addTaxonML(added_taxon, added_node, child, node,
                                  isAddedAtMidpoint,
                                  target_node2, target_dad2,
                                  len_child, len_node, len_dad);
        if (score > best_score) {
            best_score         = score;
            target_node        = target_node2;
            target_dad         = target_dad2;
            len_to_new_taxon   = len_child;
            len_to_target_node = len_node;
            len_to_target_dad  = len_dad;
        }
    }
    return best_score;
}

class BlockAllocator {
    
private:
    PhyloTree* phylo_tree;
    size_t     nptn;
    uint64_t   parsimony_block_size;
    uint64_t   lh_block_size;
    uint64_t   scale_block_size;
    int        index_parsimony;
    int        index_lh;
    
public:
    BlockAllocator(PhyloTree* tree,
                   int& parsimonyIndex, int& likelihoodIndex)
        : phylo_tree(tree), index_parsimony(parsimonyIndex), index_lh(likelihoodIndex) {
        tree->getBlockSizes( nptn, parsimony_block_size, lh_block_size, scale_block_size );
    }
    void allocateLikelihoodBlocks(double*& partial_lh, UBYTE*& scale_num) {
        partial_lh = phylo_tree->central_partial_lh + (index_lh * lh_block_size);
        scale_num  = phylo_tree->central_scale_num  + (index_lh * scale_block_size);
        ++index_lh;
    }
    void allocateParsimonyBlock(UINT*& partial_pars) {
        partial_pars = phylo_tree->central_partial_pars + (index_parsimony * parsimony_block_size);
        ++index_parsimony;
    }
    void allocateMemoryFor(PhyloNeighbor* nei) {
        if (nei->partial_lh==nullptr) {
            allocateLikelihoodBlocks(nei->partial_lh, nei->scale_num);
        }
        if (nei->partial_pars==nullptr) {
            allocateParsimonyBlock(nei->partial_pars);
        }
    }
    PhyloTree* getTree() {
        return phylo_tree;
    }
    int getLikelihoodBlockCount() const {
        return index_lh;
    }
    int getParsimonyBlockCount() const {
        return index_parsimony;
    }
    void handOverComputedState(PhyloNeighbor* from_nei, PhyloNeighbor* to_nei) {
        std::swap(to_nei->partial_lh         , from_nei->partial_lh);
        std::swap(to_nei->partial_pars       , from_nei->partial_pars);
        std::swap(to_nei->scale_num          , from_nei->scale_num);
        std::swap(to_nei->partial_lh_computed, from_nei->partial_lh_computed);
        allocateMemoryFor(from_nei);
        to_nei->partial_lh_computed = from_nei->partial_lh_computed;
        from_nei->clearComputedFlags();
    }
};

ParallelParsimonyCalculator::ParallelParsimonyCalculator(PhyloTree& phylo_tree)
    : tree(phylo_tree), task_to_start(nullptr), task_in_progress(nullptr) {}

void ParallelParsimonyCalculator::computePartialParsimony
    ( PhyloNeighbor* dad_branch, PhyloNode* dad ) {
    if (!dad_branch->isParsimonyComputed()) {
        workToDo.emplace_back(dad_branch, dad);
    }
}

int  ParallelParsimonyCalculator::computeParsimonyBranch
    ( PhyloNeighbor* dad_branch, PhyloNode* dad,
      const char* task_description ) {
    PhyloNode*     node        = dad_branch->getNode();
    PhyloNeighbor* node_branch = node->findNeighbor(dad);

    int startIndex = workToDo.size();
    computePartialParsimony(dad_branch,  dad);
    computePartialParsimony(node_branch, node);
    calculate(startIndex, task_description);
    return tree.computeParsimonyBranch(dad_branch, dad);
}

void ParallelParsimonyCalculator::calculate(int start_index,
                                            const char* task_description) {
    int stop_index = workToDo.size();
    bool tasked = (task_description!=nullptr && task_description[0]!='\0');
    if (stop_index <= start_index) {
        //Bail, if nothing to do
        return;
    }
    
    if (tasked && task_to_start==nullptr) {
        task_to_start   = task_description;
        task_in_progress = task_description;
    }
    
    //1. Find work to do at a lower level
    for (int i=stop_index-1; i>=start_index; --i) {
        auto item = workToDo[i];
        PhyloNeighbor* dad_branch = item.first;
        PhyloNode*     dad        = item.second;
        PhyloNode*     node       = dad_branch->getNode();
        FOR_EACH_PHYLO_NEIGHBOR(node, dad, it, nei) {
            computePartialParsimony(nei, node);
        }
    }
    
    //2. Do it, and then forget about it
    calculate(stop_index, "");
    workToDo.resize(stop_index);
    
    //3. Do the actual parsimony calculations at the
    //   current level (this doesn't change the content
    //   of workToDo so workToDo's contents can be
    //   processed with a parallel pointer loop).
    WorkItem* startItem = workToDo.data() + start_index;
    WorkItem* stopItem  = workToDo.data() + stop_index;
    if (task_to_start != nullptr) {
        tree.initProgress( workToDo.size(), task_to_start, "", "" );
        task_to_start = nullptr;
    }
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (WorkItem* item = startItem; item<stopItem; ++item) {
        PhyloNeighbor* dad_branch = item->first;
        PhyloNode*     dad        = item->second;
        tree.computePartialParsimony(dad_branch, dad);
    }
    if (task_in_progress != nullptr) {
        tree.trackProgress(stopItem-startItem);
    }
    workToDo.resize(start_index);
    if (tasked) {
        tree.doneProgress();
        task_in_progress = nullptr;
    }
}

class TaxonToPlace;      //A taxon being added
class TargetBranch;      //A place it could go
class TargetBranchRange; //A container of places taxa could go
class TargetBranchRef;   //A reference to a target branch via a target branch range
class PossiblePlacement; //A costing for a TaxonToPlace at a TargetBranch

class SearchHeuristic {
public:
    virtual ~SearchHeuristic() = default;
    virtual bool isPlacementWorthTrying
        ( const TaxonToPlace* taxon, const TargetBranch* target ) {
            return true;
    }
    virtual bool isGlobalSearch() {
        return true;
    }
};

class PlacementCostCalculator {
public:
    virtual void assessPlacementCost(PhyloTree& tree, const TaxonToPlace* taxon,
                                     PossiblePlacement* p) const = 0;
    virtual bool usesLikelihood() {
        return false;
    }
    virtual ~PlacementCostCalculator() = default;
};

typedef std::vector<TargetBranchRef>
    ReplacementBranchList;

class TargetBranch : public std::pair<PhyloNode*, PhyloNode*> {
    //A place where a node could be inserted, with likelihood and
    //partial parsimony determined, looking into the tree from the
    //insertion point.
    UINT*   partial_pars;
    double* partial_lh;
    UBYTE*  scale_num;
    bool    used;
    ReplacementBranchList* replacements;
    friend class TargetBranchRef;
public:
    typedef std::pair<PhyloNode*, PhyloNode*> super;
    TargetBranch(): super(nullptr, nullptr)
                    , partial_pars(nullptr), partial_lh(nullptr)
                    , scale_num(nullptr), used(false)
                    , replacements(nullptr) {}
    ~TargetBranch() {
        delete replacements;
    }
    TargetBranch(const TargetBranch& rhs): super(rhs) {
        partial_pars = rhs.partial_pars;
        partial_lh   = rhs.partial_lh;
        scale_num    = rhs.scale_num;
        used         = rhs.used;
        replacements = (rhs.replacements==nullptr) ? nullptr
        : new ReplacementBranchList(*rhs.replacements);
    }
    TargetBranch& operator= (const TargetBranch& rhs) {
        super::operator=(rhs);
        partial_pars = rhs.partial_pars;
        partial_lh   = rhs.partial_lh;
        scale_num    = rhs.scale_num;
        used         = rhs.used;
        delete         replacements;
        replacements = (rhs.replacements==nullptr) ? nullptr
        : new ReplacementBranchList(*rhs.replacements);
        return *this;
    }
    TargetBranch(BlockAllocator& allocator,
                   PhyloNode* node1, PhyloNode* node2,
                   bool likelihood_wanted)
        : super(node1, node2), replacements(nullptr) {
        allocator.allocateParsimonyBlock(partial_pars);
        if (!likelihood_wanted) {
            allocator.allocateLikelihoodBlocks(partial_lh, scale_num);
        }
        else {
            partial_lh = nullptr;
            scale_num  = nullptr;
        }
        used  = false;
    }
    void computeState(PhyloTree& phylo_tree) const {
        PhyloNeighbor* neigh1 = first->findNeighbor(second);
        PhyloNeighbor* neigh2 = second->findNeighbor(first);
        ParallelParsimonyCalculator c(phylo_tree);
        c.computePartialParsimony(neigh1, first);
        c.computePartialParsimony(neigh2, second);
        c.calculate();
        phylo_tree.computePartialParsimonyOutOfTree
            ( neigh1->partial_pars, neigh2->partial_pars, partial_pars );
        //Todo: set up local likelihood and scale_num vectors,
        //if partial_lh and scale_num are not... nullptr.
    }
    bool isUsedUp() const {
        return used;
    }
    void handOverComputedStateTo(PhyloNeighbor* nei) {
        nei->partial_pars = partial_pars;
        nei->partial_lh   = partial_lh;
        nei->scale_num    = scale_num;
        partial_pars      = nullptr;
        partial_lh        = nullptr;
        scale_num         = nullptr;
        nei->setParsimonyComputed  ( true );
        nei->setLikelihoodComputed ( nei->partial_lh != nullptr );
        used  = true;
    }
    const UINT* getParsimonyBlock() const {
        return partial_pars;
    }
    template <class T>
    void costPlacementOfTaxa(PhyloTree& tree,
                             TargetBranchRange* targets,
                             size_t targetNumber,
                             T* candidateStart,
                             T* candidateStop,
                             SearchHeuristic*   heuristic,
                             PlacementCostCalculator* calculator,
                             bool isFirstTargetBranch) const;
    void takeOwnershipOfReplacementVector(ReplacementBranchList* branches) {
        replacements = branches;
    }
    ReplacementBranchList* getReplacements() {
        return replacements;
    }
};

class TargetBranchRef;

class TargetBranchRange : public vector<TargetBranch> {
public:
    typedef  vector<TargetBranch> super;
    TargetBranchRange(PhyloTree& phylo_tree, BlockAllocator& b,
                   PlacementCostCalculator& calculator): super() {
        PhyloNodeVector v1, v2;
        phylo_tree.getBranches(v1, v2);
        reserve(v1.size());
        if ( verbose_mode >= VB_DEBUG ) {
            std::stringstream s1;
            s1 << "TargetBranchRange will have " << v1.size() << " entries";
            phylo_tree.logLine(s1.str());
        }
        for (int i=0; i<v1.size(); ++i) {
            emplace_back(b, v1[i], v2[i], calculator.usesLikelihood());
        }
    }
    void removeUsed() {
        int w=0;
        for (int r=0; r<size(); ++r) {
            if (!at(r).isUsedUp()) {
                at(w) = at(r);
                ++w;
            }
        }
        resize(w);
    }
    TargetBranchRef addNewRef(BlockAllocator& allocator,
                              PhyloNode* node1, PhyloNode* node2,
                              bool likelihood_wanted);
};

class TargetBranchRef {
private:
    TargetBranchRange* target_range;
    size_t             target_index;
public:
    TargetBranchRef(): target_range(nullptr), target_index(0) {}
    TargetBranchRef(const TargetBranchRef& r)
        : target_range(r.target_range), target_index(r.target_index) {}
    TargetBranchRef(TargetBranchRange* range, size_t index)
        : target_range(range), target_index(index) {}
    TargetBranchRef& operator=(const TargetBranchRef& r) = default;
    bool isUsedUp() const {
        if (target_range == nullptr) {
            return true;
        }
        return target_range->at(target_index).isUsedUp();
    }
    PhyloNode* getFirst() const {
        if (target_range == nullptr) {
            return nullptr;
        }
        return target_range->at(target_index).first;
    }
    PhyloNode* getSecond() const {
        if (target_range == nullptr) {
            return nullptr;
        }
        return target_range->at(target_index).second;
    }
    TargetBranch* getTarget() const {
        if (target_range == nullptr) {
            return nullptr;
        }
        return &target_range->at(target_index);
    }
    size_t getTargetIndex() const {
        return target_index;
    }
};

TargetBranchRef TargetBranchRange::addNewRef(BlockAllocator& allocator,
                                             PhyloNode* node1, PhyloNode* node2,
                                             bool likelihood_wanted) {
    emplace_back(allocator, node1, node2, likelihood_wanted);
    back().computeState(*allocator.getTree());
    return TargetBranchRef(this, size()-1);
}

class PossiblePlacement {
public:
    TargetBranchRef  target_branch;
    const PhyloNode* node1;          //used to check if the insertion point
    const PhyloNode* node2;          //still exists (when about to insert).
    double           score;          //score (the likelihood, or minus
                                     //the parsimony score)
    double           lenToNewTaxon;  //(best scoring) length of the edge between
                                     //new_taxon and added_node
    double           lenToNode1;     //(best scoring) length of edge between
                                     //target_dad and added_node
    double           lenToNode2;     //(best scoring) length of edge between
                                     //target_child and added_node
    TargetBranchRef  replacement_start;
    TargetBranchRef  replacement_stop;
    
    PossiblePlacement(): node1(nullptr), node2(nullptr)
                        , score(0), lenToNewTaxon(-1)
                        , lenToNode1(0), lenToNode2(0) {
    }
    PossiblePlacement& operator= (const PossiblePlacement & rhs) = default;
    bool operator < (const PossiblePlacement& rhs) const {
        return score < rhs.score;
    }
    bool operator <= (const PossiblePlacement& rhs) const {
        return score <= rhs.score;
    }
    void setTargetBranch(TargetBranchRange* targetRange, size_t index) {
        target_branch = TargetBranchRef(targetRange, index);
        node1 = target_branch.getFirst();
        node2 = target_branch.getSecond();
    }
    void setTargetBranch(TargetBranchRef& branch_ref) {
        target_branch = branch_ref;
        node1 = target_branch.getFirst();
        node2 = target_branch.getSecond();
    }
    bool canStillUse() const {
        return !target_branch.isUsedUp();
    }
    TargetBranch* getTarget() const {
        return target_branch.getTarget();
    }
    size_t getTargetIndex() const {
        return target_branch.getTargetIndex();
    }
    void forget() {
        target_branch = TargetBranchRef(nullptr, 0);
    }
};

class TaxonToPlace {
    //A taxon that could be added to a tree.
protected:
    PossiblePlacement bestPlacement;
public:
    int               taxonId;
    std::string       taxonName;
    bool              inserted;      //true if this taxon has been inserted
    PhyloNode*        new_leaf;      //leaf
    PhyloNode*        new_interior;  //interior
    const UINT*       partial_pars;  //partial parsimony for new leaf, seen
                                     //from the new interior
    
    TaxonToPlace() : bestPlacement(), taxonId(-1), taxonName()
                   , inserted(false), new_leaf(nullptr), new_interior(nullptr)
                   , partial_pars(nullptr) {
    }

    TaxonToPlace(const TaxonToPlace& rhs) = default; //copy everything!
    
    TaxonToPlace(BlockAllocator& ba, int id, std::string name)
        : taxonId(id), taxonName(name), inserted(false) {
        PhyloTree* phylo_tree = ba.getTree();
        new_leaf     = phylo_tree->newNode(taxonId, taxonName.c_str());
        new_interior = phylo_tree->newNode();
        new_interior->addNeighbor(new_leaf, -1 );
        new_leaf->addNeighbor(new_interior, -1 );
        PhyloNeighbor* nei = new_interior->firstNeighbor();
        ba.allocateMemoryFor(nei);
        phylo_tree->computePartialParsimony(nei, new_interior);
        partial_pars = new_interior->firstNeighbor()->partial_pars;
    }
    
    virtual ~TaxonToPlace() = default;
        
    const UINT* getParsimonyBlock() const {
        return partial_pars;
    }
    
    int getTaxonId() const {
        return taxonId;
    }
    
    virtual size_t considerPlacements(const PossiblePlacement* placements, size_t count) {
        size_t bestI = 0;
        for (size_t i=1; i<count; ++i) {
            if (placements[bestI].score < placements[i].score ) {
                bestI = i;
            }
        }
        bestPlacement = placements[bestI];
        return bestI;
    }
    
    virtual bool considerAdditionalPlacement(const PossiblePlacement& placement) {
        bool best = (!bestPlacement.canStillUse() || placement.score < bestPlacement.score);
        if (best) {
            bestPlacement = placement;
        }
        return best;
    }
    
    virtual const PossiblePlacement& getBestPlacement() const {
        return bestPlacement;
    }

    virtual bool canInsert() const {
        return bestPlacement.canStillUse();
    }

    void findPlacement(PhyloTree& phylo_tree,
                       TargetBranchRange* range,
                       size_t insertion_point_count,
                       SearchHeuristic* heuristic,
                       const PlacementCostCalculator* calculator) {
        //
        //Note: for now, heuristic is ignored...
        //
        std::vector<PossiblePlacement> placements;
        placements.resize(insertion_point_count);
        PossiblePlacement* placementArray = placements.data();
        
        if ( verbose_mode >= VB_DEBUG ) {
            std::stringstream s1;
            s1 << "Scoring " << this->taxonName;
            phylo_tree.logLine(s1.str());
        }

        auto rangeStart = range->data();
        #ifdef _OPENMP
        #pragma omp parallel for
        #endif
        for (int i=0; i<insertion_point_count; ++i) {
            if (heuristic->isPlacementWorthTrying(this, rangeStart+i)) {
                PossiblePlacement* p = placementArray + i;
                p->setTargetBranch(range, i);
                calculator->assessPlacementCost(phylo_tree, this, p);
            }
        }
        phylo_tree.trackProgress(insertion_point_count);
        
        auto bestI = considerPlacements(placements.data(), placements.size());
        
        if ( verbose_mode >= VB_MED ) {
            std::stringstream s3;
            s3 << "Best (lowest) score for " << this->taxonName << " was " << bestPlacement.score << " at place " << bestI;
            phylo_tree.logLine(s3.str());
        }
        inserted      = false;
    }
    bool operator < (const TaxonToPlace& rhs) const {
        return bestPlacement.score < rhs.bestPlacement.score; //lowest score is best
    }
    bool operator <= (const TaxonToPlace& rhs) const {
        return bestPlacement.score <= rhs.bestPlacement.score; //lowest score is best
    }
    void insertIntoTree(PhyloTree& phylo_tree, BlockAllocator& b,
                        TargetBranchRange& dest,
                        PlacementCostCalculator& calculator) {
        //
        //Assumes, canInsert() returned true, and the tree has not
        //been modified in the meantime.
        //
        PhyloNode* node_1    = const_cast<PhyloNode*>(bestPlacement.node1);
        PhyloNode* node_2    = const_cast<PhyloNode*>(bestPlacement.node2);
        TargetBranch* target = bestPlacement.target_branch.getTarget();
        
        new_interior->findNeighbor(new_leaf)->length = bestPlacement.lenToNewTaxon;
        new_leaf->findNeighbor(new_interior)->length = bestPlacement.lenToNewTaxon;
        new_interior->addNeighbor(node_1,   bestPlacement.lenToNode1 );
        new_interior->addNeighbor(node_2,   bestPlacement.lenToNode2 );
        
        b.handOverComputedState( node_1->findNeighbor(node_2), new_interior->findNeighbor(node_2) );
        b.handOverComputedState( node_2->findNeighbor(node_1), new_interior->findNeighbor(node_1) );
        target->handOverComputedStateTo( new_leaf->findNeighbor(new_interior));
        
        node_1->updateNeighbor (node_2, new_interior, bestPlacement.lenToNode1);
        node_2->updateNeighbor (node_1, new_interior, bestPlacement.lenToNode2);
        
        //Todo: Are these two really needed?  Won't they be computed later anyway?
        //As and when needed?
        phylo_tree.computeParsimonyBranch(node_1->findNeighbor(new_interior), node_1);
        phylo_tree.computeParsimonyBranch(node_2->findNeighbor(new_interior), node_2);
        
        //Todo: redo likelihood too
        inserted = true;
        bool  likelihood_needed = calculator.usesLikelihood();
        ReplacementBranchList* reps = new ReplacementBranchList;
        reps->emplace_back ( dest.addNewRef(b, new_interior, node_1  , likelihood_needed ) );
        reps->emplace_back ( dest.addNewRef(b, new_interior, node_2  , likelihood_needed ) );
        reps->emplace_back ( dest.addNewRef(b, new_interior, new_leaf, likelihood_needed ) );
        bestPlacement.target_branch.getTarget()->takeOwnershipOfReplacementVector(reps);
    }
    virtual void forgetGazumpedPlacements() {
        bestPlacement.forget();
    }
    bool insertNearby(PhyloTree& phylo_tree, BlockAllocator& b,
                      TargetBranchRange& dest,
                      PlacementCostCalculator& calculator ) {
        TargetBranch* blocked_target = bestPlacement.getTarget();
        forgetGazumpedPlacements();
        std::vector<PossiblePlacement> placements;
        assessNewTargetBranches(phylo_tree, calculator,
                                blocked_target, placements);
        for (auto it = placements.begin(); it != placements.end(); ++it) {
            considerAdditionalPlacement(*it);
        }
        if (!canInsert()) {
            return false;
        }
        insertIntoTree(phylo_tree, b, dest, calculator);
        return true;
    }
    void assessNewTargetBranches(PhyloTree& phylo_tree, PlacementCostCalculator& calculator,
                                 TargetBranch* tb, std::vector<PossiblePlacement>& scores) {
        if (tb==nullptr  ) {
            return;
        }
        auto reps = tb->getReplacements();
        if (reps==nullptr) {
            return;
        }
        std::vector< ReplacementBranchList* > stack;
        stack.push_back(reps);
        while (!stack.empty()) {
            reps = stack.back();
            stack.pop_back();
            for (auto it=reps->begin(); it!=reps->end(); ++it) {
                if (it->isUsedUp()) {
                    stack.push_back(it->getTarget()->getReplacements());
                } else {
                    PossiblePlacement p;
                    p.setTargetBranch(*it);
                    calculator.assessPlacementCost(phylo_tree, this, &p);
                    scores.emplace_back(p);
                }
            }
        }
    }
};

template <class T=TaxonToPlace> class TaxaToPlace: public std::vector<T>
{
public:
    typedef std::vector<T> super;
    explicit TaxaToPlace(size_t reservation) {
        super::reserve(reservation);
    }
    virtual ~TaxaToPlace() = default;
};

template <class T> //T=TaxonToPlace, usually
void TargetBranch::costPlacementOfTaxa
    (PhyloTree&         phylo_tree,
     TargetBranchRange* targets,
     size_t             targetNumber,
     T*                 candidateStart,
     T*                 candidateStop,
     SearchHeuristic*   heuristic,
     PlacementCostCalculator* calculator,
     bool               isFirstTargetBranch) const {
    double candidateCount = candidateStop - candidateStart;
    #ifdef _OPENMP
    #pragma omp parallel for
    #endif
    for (T* candidate = candidateStart;
         candidate < candidateStop; ++candidate) {
        if (heuristic->isPlacementWorthTrying(candidate, this)) {
            PossiblePlacement p;
            p.setTargetBranch(targets, targetNumber);
            calculator->assessPlacementCost(phylo_tree, candidate, &p);
            candidate->considerAdditionalPlacement(p);
        }
    }
    phylo_tree.trackProgress(candidateCount);
}

class ParsimonyCostCalculator : public PlacementCostCalculator {
    virtual void assessPlacementCost(PhyloTree& phylo_tree,
                                     const TaxonToPlace* taxon,
                                     PossiblePlacement* placement)  const {
        auto target = placement->getTarget();
        int score;
        phylo_tree.computeParsimonyOutOfTree( target->getParsimonyBlock(),
                                              taxon->getParsimonyBlock(),
                                              &score);
        placement->score = score;
        if ( verbose_mode >= VB_MAX ) {
            std::stringstream s2;
            s2  << "Parsimony score for taxon " << taxon->getTaxonId() << " at "
                << placement->getTargetIndex() << " would be " << score;
            phylo_tree.logLine(s2.str());
        }
    }
};

class TaxonCleaner {
public:
    TaxonCleaner() = default;
    virtual ~TaxonCleaner() = default;
    void cleanUpAfterTaxonPlacement(const TaxonToPlace& taxon, PhyloTree* tree) {
    }
};
class BatchCleaner {
public:
    BatchCleaner() = default;
    virtual ~BatchCleaner() = default;
    template <class T> void cleanUpAfterBatch(TaxaToPlace<T>& taxa,
                           int firstTaxon, int lastTaxon,
                           PhyloTree* tree) {
        if (VB_MIN <= verbose_mode) {
            std::stringstream s;
            s << "Processed batch of "
              << (lastTaxon - firstTaxon) << " taxa";
            tree->logLine(s.str() );
        }
    }
};

class GlobalCleaner {
public:
    GlobalCleaner() = default;
    virtual ~GlobalCleaner() = default;
    void cleanUpAfterPlacement(PhyloTree* tree) {
    }
};

class LessFussyTaxon: public TaxonToPlace {
protected:
    std::vector<PossiblePlacement> placementStore;
    static const size_t max_placements_to_keep = 5;
    //Note: this had better not be more than about 10, because
    //      if it were large, you'd want to maintain a heap.
public:
    typedef TaxonToPlace super;
    LessFussyTaxon() : super() { }
    LessFussyTaxon(const LessFussyTaxon& rhs) = default; //copy everything!
    LessFussyTaxon(BlockAllocator& ba, int id, std::string name)
        : super(ba, id, name) {
    }
    virtual size_t considerPlacements(const PossiblePlacement* placements, size_t count) {
        placementStore.clear();
        size_t            bestI = 0;
        for (int i=0; i < count; ++i) {
            considerAdditionalPlacement(placements[i]);
            if (i==0 || placements[i] < placements[bestI]) {
                bestI    = i;
            }
        }
        return bestI;
    }
    virtual bool considerAdditionalPlacement(const PossiblePlacement& placement) {
        size_t placement_count = placementStore.size();
        if (max_placements_to_keep <= placement_count) {
            if (placementStore.back() < placement) {
                return false;
            }
            else {
                placementStore.pop_back();
            }
        }
        placementStore.emplace_back(PossiblePlacement());
        size_t sweep = placementStore.size() - 1;
        while ( 0 < sweep && placement < placementStore[sweep-1]) {
            placementStore[sweep] = placementStore[sweep - 1];
            --sweep;
        }
        placementStore[sweep] = placement;
        if (sweep == 0) {
            bestPlacement = placement;
            return true;
        }
        return false;
    }
    virtual void forgetGazumpedPlacements() {
        size_t w = 0;
        for (size_t r = 0; r < placementStore.size(); ++r) {
            if (placementStore[r].canStillUse()) {
                if (w < r) {
                    placementStore[w] = placementStore[r];
                }
                ++w;
            }
        }
        placementStore.resize(w);
        if (w == 0) {
            bestPlacement.forget();
        }
        else {
            bestPlacement = placementStore[0];
        }
    }
};

namespace {
    TaxonCleaner* getTaxonCleaner() {
        auto localCleanup = getLocalCleanupAlgorithm();
        switch (localCleanup) {
            default:
                break;
        }
        return new TaxonCleaner();
    }
    BatchCleaner* getBatchCleaner() {
        auto batchCleanup = getBatchCleanupAlgorithm();
        switch (batchCleanup) {
            default:
                break;
        }
        return new BatchCleaner();
    }
    GlobalCleaner* getGlobalCleaner() {
        auto globalCleanup = getGlobalCleanupAlgorithm();
        switch (globalCleanup) {
            default:
                break;
        }
        return new GlobalCleaner();
    }
}

void PhyloTree::removeSampleTaxaIfRequested() {
    size_t nseq = aln->getNSeq();
    size_t countOfTaxaToRemove = getNumberOfTaxaToRemove(nseq);
    if (0<countOfTaxaToRemove) {
        map<string, Node*> mapNameToNode;
        getMapOfTaxonNameToNode(nullptr, nullptr, mapNameToNode);
        size_t r = 0;
        for (size_t seq = 0; seq < nseq; ++seq) {
            r += countOfTaxaToRemove;
            if ( r >= nseq ) {
                r -= nseq;
                string seq_name = aln->getSeqName(seq);
                auto it = mapNameToNode.find(seq_name);
                if (it!=mapNameToNode.end()) {
                    Node* node = it->second;
                    auto newName = node->name + "_Removed";
                    if (mapNameToNode.find(newName) == mapNameToNode.end()) {
                        node->name = newName;
                        mapNameToNode[newName] = node;
                    }
                }
            }
        }
    }
}

double PhyloTree::taxaAdditionWorkEstimate(size_t newTaxaCount, size_t taxaPerBatch,
                                           size_t insertsPerBatch) {
    if ( newTaxaCount <= taxaPerBatch || taxaPerBatch == 0 ) {
        if ( newTaxaCount <= insertsPerBatch || insertsPerBatch == 0 ) {
            return 3.0 * newTaxaCount * leafNum;
        }
        return 3.0 * newTaxaCount * leafNum * newTaxaCount / insertsPerBatch;
    }
    size_t batchesThisPass  = newTaxaCount / taxaPerBatch;
    double workThisPass     = batchesThisPass * taxaPerBatch * leafNum;
    double progressThisPass = batchesThisPass * insertsPerBatch;
    //Optimistic if inserts = 100% and batches are large.
    return (3.0 * workThisPass / progressThisPass) * newTaxaCount;
}

bool PhyloTree::shouldPlacementUseSankoffParsimony() const {
    return getCostFunction() == SANKOFF_PARSIMONY;
}

bool PhyloTree::shouldPlacementUseLikelihood() const {
    CostFunction       costFunction    = getCostFunction();
    return ( costFunction != MAXIMUM_PARSIMONY
             && costFunction != SANKOFF_PARSIMONY);
}

void logInsert(PhyloTree* tree, Params& params, CostFunction costFunction,
               size_t totalInsertCount, const char* verb,
               TaxonToPlace & c, const char * where) {
    if (( verbose_mode >= VB_MIN && !params.suppress_list_of_sequences)
        || verbose_mode >= VB_MED ) {
        stringstream s;
        s << totalInsertCount << ". " << verb << " "
            << c.taxonName << " " << where << ". It had ";
        const PossiblePlacement& p = c.getBestPlacement();
        if (costFunction==MAXIMUM_PARSIMONY ||
            costFunction==SANKOFF_PARSIMONY) {
            s << "parsimony score " << (int)(p.score);
        } else {
            s << "likelihood score " << p.score;
        }
        s << " (and path lengths " << p.lenToNode1
            << ", " << p.lenToNode2 << ", " << p.lenToNewTaxon << ")";
        tree->logLine(s.str());
    }
}

SearchHeuristic* getSearchHeuristic() {
    return new SearchHeuristic;
}

typedef TaxonToPlace TaxonTypeInUse;
//typedef LessFussyTaxon TaxonTypeInUse;

#define NEW_TAXON_MAJOR (0)

struct PlacementRun {
public:
    CostFunction             costFunction; //parsimony or likelihood
    SearchHeuristic*         heuristic; //global? or localized?
    TaxonCleaner*            taxon_cleaner;
    BatchCleaner*            batch_cleaner;
    GlobalCleaner*           global_cleaner;
    PlacementCostCalculator* calculator;

    PlacementRun(): costFunction(getCostFunction())
        , heuristic(getSearchHeuristic())
        , taxon_cleaner(getTaxonCleaner())
        , batch_cleaner(getBatchCleaner())
        , global_cleaner(getGlobalCleaner())
        , calculator(new ParsimonyCostCalculator()) {
       costFunction = getCostFunction();
    }
    ~PlacementRun() {
        delete global_cleaner;
        delete batch_cleaner;
        delete taxon_cleaner;
        delete heuristic;
        delete calculator;
    }
};

void PhyloTree::addNewTaxaToTree(const IntVector& taxaIdsToAdd) {
    //
    //Assumes: The tree is rooted.
    //
    Params&            params          = Params::getInstance();
    size_t             taxaPerBatch    = getTaxaPerBatch(taxaIdsToAdd.size());
                                         //Must be 1 or more
    size_t             insertsPerBatch = getInsertsPerBatch(taxaIdsToAdd.size(), taxaPerBatch);
                                         //Must be 1 or more
    PlacementRun pr;
    deleteAllPartialLh();
    
    if ( taxaPerBatch == 1 && pr.heuristic->isGlobalSearch()  &&
        ( pr.costFunction == MAXIMUM_PARSIMONY || pr.costFunction == SANKOFF_PARSIMONY ) ) {
        //For now, we might as well use the existing step-wise
        //parsimony stuff for adding to a constraint tree, eh?
        //Since, for now, it is a lot faster.
        constraintTree.readConstraint(*this);
        //clearing all the nodes...
        freeNode();
        root = nullptr;
        cout << "Creating fast initial parsimony tree by random order stepwise addition..." << endl;
        double start = getRealTime();
        double score = computeParsimonyTree(params.out_prefix, aln, randstream);
        cout << getRealTime() - start << " seconds, parsimony score: " << score
            << " (based on " << aln->num_parsimony_sites << " sites)"<< endl;
        
        //Note that this score tends to disagree.
        double parsimonyStart = getRealTime();
        clearAllPartialParsimony(false);
        double parsimonyScore = computeParsimony("Recalculating parsimony score");
        LOG_LINE( VB_MED, "Recalculated parsimony score " << parsimonyScore
                    << " (recalculation cost " << (getRealTime() - parsimonyStart) << " sec)" );

        finishUpAfterTaxaAddition();
        return;
    }
    bool trackLikelihood = shouldPlacementUseLikelihood();
    
    size_t   extra_parsimony_blocks = leafNum * 2 - 4;
    size_t   extra_lh_blocks        = trackLikelihood
                                    ? (leafNum * 4 - 4 + taxaIdsToAdd.size())
                                    : 0;
    int      index_parsimony        = 0;
    int      index_lh               = 0;
    ensurePartialLHIsAllocated(extra_parsimony_blocks, extra_lh_blocks);
    initializeAllPartialLh(index_parsimony, index_lh, trackLikelihood);
    if (pr.costFunction == SANKOFF_PARSIMONY) {
        computeTipPartialParsimony();
    }
    BlockAllocator allocator(this, index_parsimony, index_lh);
    
    LOG_LINE ( VB_MED, "After overallocating lh blocks, index_lh was " << allocator.getLikelihoodBlockCount() );
    if (VB_MED <= verbose_mode) {
        curScore = computeLikelihood();
        LOG_LINE ( VB_MED, "Likelihood score before insertions was " << curScore );
        #if (0)
            curScore = optimizeAllBranches(2);
            LOG_LINE ( VB_MED, "Optimized likelihood score before insertions was " << curScore);
        #endif
    }
    LOG_LINE ( VB_MED, "Batch size is " << taxaPerBatch
              << " and the number of inserts per batch is " << insertsPerBatch);
    
    double setUpStartTime = getRealTime();
    size_t newTaxaCount = taxaIdsToAdd.size();
    
    TaxaToPlace<TaxonTypeInUse> candidates(newTaxaCount);
    LOG_LINE ( VB_DEBUG, "Before allocating TaxonToPlace array"
              << ", index_lh was " << allocator.getLikelihoodBlockCount() );
    for (size_t i=0; i<newTaxaCount; ++i) {
        int         taxonId   = taxaIdsToAdd[i];
        std::string taxonName = aln->getSeqName(taxonId);
        candidates.emplace_back(allocator, taxonId, taxonName);
    }
    LOG_LINE ( VB_DEBUG, "After allocating TaxonToPlace"
              << ", index_lh was " << allocator.getLikelihoodBlockCount()
              << ", index_pars was " << allocator.getParsimonyBlockCount());

    TargetBranchRange targets(*this, allocator, *pr.calculator);
    LOG_LINE ( VB_DEBUG, "After allocating TargetBranchRange"
               << ", index_lh was " << allocator.getLikelihoodBlockCount()
               << ", index_pars was " << allocator.getParsimonyBlockCount());
    LOG_LINE ( VB_MIN, "Set up time was " << (getRealTime() - setUpStartTime) << " sec");
    
    
    double estimate = taxaAdditionWorkEstimate
                      ( newTaxaCount, taxaPerBatch, insertsPerBatch );
    size_t totalInsertCount     = 0;
    size_t blockedInsertCount   = 0;
    double timeSpentOnRefreshes = 0.0; //Time spent recalculating parsimony &/or likelihood
                                       //for the entire tree
    double timeSpentOnSearches  = 0.0;
    double timeSpentOnInserts   = 0.0;
    initProgress(estimate, "Adding new taxa to tree", "", "");
    while (0<newTaxaCount) {
        if (newTaxaCount<taxaPerBatch) {
            taxaPerBatch = newTaxaCount;
        }
        size_t batchStart=0;
        for (; batchStart+taxaPerBatch <= newTaxaCount; batchStart+=taxaPerBatch) {
            timeSpentOnRefreshes -= getRealTime();
            if (trackLikelihood) {
                clearAllPartialLH(false);
                clearAllScaleNum(false);
                double likelihoodScore = computeLikelihood();
                LOG_LINE( VB_MIN, "Log-likelihood is currently " << likelihoodScore);
            }
            size_t batchStop = batchStart + taxaPerBatch;
            TargetBranch* pointStart = targets.data();
            TargetBranch* pointStop  = pointStart + targets.size();
            clearAllPartialParsimony(false);
            for (TargetBranch* point = pointStart; point<pointStop; ++point) {
                point->computeState(*this);
            }
            timeSpentOnRefreshes += getRealTime();
            timeSpentOnSearches -= getRealTime();
            TaxonTypeInUse* candidateStart = candidates.data() + batchStart;
            TaxonTypeInUse* candidateStop  = candidates.data() + batchStop;
#if (NEW_TAXON_MAJOR)
            for (TaxonTypeInUse* c = candidateStart; c<candidateStop; ++c) {
                LOG_LINE(VB_DEBUG, "Scoring ... " << c.taxonName);
                c->findPlacement(*this, targets.data(), targets.size(),
                                heuristic, calculator);
                PossiblePlacement& p = c.bestPlacement;
                LOG_LINE(VB_DEBUG, "Scored " << p.score << " for placement"
                         << " of " << c.taxonName << " with lengths "
                         << p.lenToNode1 << ", " << p.lenToNode2 << ", " << p.lenToNewTaxon);
            }
#else //INSERTION_POINT_MAJOR
            for (TargetBranch* point = pointStart; point<pointStop; ++point) {
                point->costPlacementOfTaxa(*this, &targets, point-pointStart,
                                           candidateStart, candidateStop,
                                           pr.heuristic, pr.calculator,
                                           point==pointStart);
            }
#endif
            timeSpentOnSearches += getRealTime();
            insertsPerBatch      = getInsertsPerBatch(taxaIdsToAdd.size(), batchStop-batchStart);
            size_t insertStop    = batchStart + insertsPerBatch;
            std::sort( candidates.begin() + batchStart, candidates.begin() + batchStop);
            if (batchStop <= insertStop) {
                insertStop = batchStop; //Want them all
            }
            timeSpentOnInserts -= getRealTime();
            size_t insertCount = 0;
            for ( size_t i = batchStart; i<insertStop; ++i) {
                TaxonToPlace& c = candidates[i];
                if (c.canInsert()) {
                    ++insertCount;
                    ++totalInsertCount;
                    c.insertIntoTree(*this, allocator, targets, *pr.calculator);
                    logInsert(this, params, pr.costFunction, totalInsertCount,
                              "Inserted", c, "at its preferred branch");
                } else {
                    //Another candidate taxon has gotten there first
                    ++blockedInsertCount;
                    ++insertCount;
                    ++totalInsertCount;
                    c.insertNearby(*this, allocator, targets, *pr.calculator);
                    logInsert(this, params, pr.costFunction, totalInsertCount,
                              "Inserted", c, "near its preferred branch");
                }
                pr.taxon_cleaner->cleanUpAfterTaxonPlacement(c, this);
            }
            timeSpentOnInserts += getRealTime();
            if ( 1 < batchStop - batchStart ) {
                LOG_LINE ( VB_MED,  "Inserted " << (insertCount)
                          << " out of a batch of " << (batchStop - batchStart) << "." );
            }
            pr.batch_cleaner->cleanUpAfterBatch(candidates, batchStart, batchStop, this);
            if (trackLikelihood) {
                fixNegativeBranch();
            }
            if (insertCount == 0) {
                outError("No taxa inserted in batch");
                break;
            }
        } //batches of items
        
        targets.removeUsed();
        //Remove all the candidates that we were able to place
        std::vector<TaxonTypeInUse> oldCandidates;
        std::swap(oldCandidates, candidates);
        //1. Any candidates not considered this time go to the
        //   first batch to consider in the next pass.
        for (size_t r=batchStart; r<newTaxaCount; ++r) {
            candidates.emplace_back(oldCandidates[r]);
        }
        //2. Any candidates that were considered, but were not
        //   inserted, are to be considered in the next pass.
        for (size_t r=0; r<batchStart; ++r) {
            if (!oldCandidates[r].inserted) {
                //Keep this one to be considered next time
                candidates.emplace_back(oldCandidates[r]);
            }
        }
        newTaxaCount = candidates.size();
        insertsPerBatch = getInsertsPerBatch(taxaIdsToAdd.size(), taxaPerBatch);
        auto workLeft   = taxaAdditionWorkEstimate
                          ( newTaxaCount, taxaPerBatch, insertsPerBatch );
        this->progress->setWorkRemaining(workLeft);
        LOG_LINE ( VB_MAX, "At the end of this pass, index_lhs was "
                  << allocator.getLikelihoodBlockCount() << ", index_pars was "
                  << allocator.getParsimonyBlockCount());
    }
    doneProgress();
    
    LOG_LINE ( VB_MED, "Tidying up tree after inserting taxa.");
    pr.global_cleaner->cleanUpAfterPlacement(this);
    
    LOG_LINE ( VB_MIN, "Time spent on refreshes was "      << timeSpentOnRefreshes << " sec");
    LOG_LINE ( VB_MIN, "Time spent on searches was "       << timeSpentOnSearches << " sec");
    LOG_LINE ( VB_MIN, "Time spent on actual inserts was " << timeSpentOnInserts << " sec");
    LOG_LINE ( VB_MIN, "Total number of blocked inserts was " << blockedInsertCount );
    LOG_LINE ( VB_MED, "At the end of addNewTaxaToTree, index_lhs was "
              << allocator.getLikelihoodBlockCount() << ", index_pars was "
              << allocator.getParsimonyBlockCount() << ".");
    if (!trackLikelihood) {
        fixNegativeBranch();
    }
    finishUpAfterTaxaAddition();
}

void PhyloTree::finishUpAfterTaxaAddition() {
    initializeTree();
    deleteAllPartialLh();
    initializeAllPartialLh();
    LOG_LINE ( VB_MED, "Number of leaves " << this->leafNum
              << ", of nodes " << this->nodeNum );
    CostFunction costFunction = getCostFunction();
    if (costFunction==MAXIMUM_LIKELIHOOD_ANYWHERE ||
        costFunction==MAXIMUM_LIKELIHOOD_MIDPOINT) {
        double score = optimizeAllBranches();
        LOG_LINE ( VB_MIN, "After optimizing, likelihood score was " << score );
    }
}
