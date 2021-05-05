/*
 * mp2.cpp
 *
 *  Created on: May 3, 2021
 *      Author: Diep Thi Hoang (diep.thi.hoang@gmail.com)
 */


#include "phylotree.h"
#include "iqtree.h"
#include "mp2.h"
#include "alignment/alignment.h"

// Trivial score test: ./iqtree2-mpi -s example.phy -t example.phy.treefile --spr-unit-test 1 -redo
void doSPRUnitTest(Params &params) {
    Alignment alignment(params.aln_file, params.sequence_type, params.intype, params.model_name);

    // Diep: 2021-05-04
    // The following line is NEW and EXTREMELY IMPORTANT for iqtree2 parsimony.
    // I'm wondering ... Does everyone know?
    alignment.orderPatternByNumChars(PAT_VARIANT);

    IQTree * ptree;
    ptree = new IQTree(&alignment);
    ptree->readTree(params.user_file.c_str(), params.is_rooted); // Read user tree
    ptree->setAlignment(&alignment); // IMPORTANT: Always call setAlignment() after readTree()
	ptree->drawTree(cout);

//    if (ptree->rooted)
//        ptree->convertToUnrooted();
//    ptree->initCostMatrix(CM_UNIFORM);
//    ptree->setParsimonyKernel(params.SSE);

	ptree->initializeAllPartialPars();
	cout << "Parsimony score by iqtree2 is: "
			<< ptree->computeParsimony()
            << endl;
	delete ptree;
}


