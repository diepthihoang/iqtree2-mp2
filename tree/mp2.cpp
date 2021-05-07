/*
 * mp2.cpp
 *
 *  Created on: May 3, 2021
 *      Author: Diep Thi Hoang (diep.thi.hoang@gmail.com)
 */


#include "mp2.h"
#include "phylotree.h"
#include "iqtree.h"
#include "alignment/alignment.h"


// Run unit test for parsimony score:
// ./iqtree2-mpi -s example.phy -t example.phy.treefile --test-pars-score 1 -redo
// --test-pars-score 1 <1|2> where 1 is for PLL core, 2 is for IQTREE2 core
void ParsimonyTest::runUnitTestParsimonyScore(Params &params) {
    if(params.unit_test_pars_score == 2){
        cout << "Test using iqtree2/dev parsimony computation:" << endl;
        Alignment alignment(params.aln_file, params.sequence_type, params.intype, params.model_name);

        // Diep: 2021-05-04
        // The following line is NEW and EXTREMELY IMPORTANT for iqtree2 parsimony.
        // I'm wondering ... Does everyone know?
        alignment.orderPatternByNumChars(PAT_VARIANT);

        IQTree * ptree;
        ptree = new IQTree(&alignment);
        ptree->readTree(params.user_file.c_str(), params.is_rooted); // Read user tree
        ptree->setAlignment(&alignment); // IMPORTANT: Always call setAlignment() after readTree()
        ptree->setParams(&params);
        // ptree->drawTree(cout);

    //    if (ptree->rooted)
    //        ptree->convertToUnrooted();
    //    ptree->initCostMatrix(CM_UNIFORM);
    //    ptree->setParsimonyKernel(params.SSE);

        ptree->initializeAllPartialPars();
        cout << "Parsimony score (by iqtree2 kernel) is: "
                << ptree->computeParsimony()
                << endl;
        delete ptree;
    }else if(params.unit_test_pars_score == 1){
        cout << "Test using PLL parsimony computation:" << endl;
        Alignment alignment(params.aln_file, params.sequence_type, params.intype, params.model_name);
        // alignment.orderPatternByNumChars(PAT_VARIANT);
        IQTree * ptree;
        ptree = new IQTree(&alignment);
        ptree->readTree(params.user_file.c_str(), params.is_rooted); // Read user tree
        ptree->setAlignment(&alignment); // IMPORTANT: Always call setAlignment() after readTree()
        ptree->setParams(&params);
        
        ptree->initializePLL(params); // Diep: This one modifies topo?
        allocateParsimonyDataStructures(ptree->pllInst, ptree->pllPartitions);

    	// string tree_string = ptree->getTreeString();
        stringstream tree_stream;
        ptree->printTree(tree_stream); // ptree->printTree(tree_stream, WT_SORT_TAXA);
        string tree_string = tree_stream.str();
        cout << "\niqtree2 getTreeString = " << tree_string << endl;
    
        // tree_string = "(LngfishAu:-1.000000000000000,(LngfishSA:-1.000000000000000,LngfishAf:-1.000000000000000):-1.000000000000000,(Frog:-1.000000000000000,((Turtle:-1.000000000000000,((Sphenodon:-1.000000000000000,Lizard:-1.000000000000000):-1.000000000000000,(Crocodile:-1.000000000000000,Bird:-1.000000000000000):-1.000000000000000):-1.000000000000000):-1.000000000000000,((((Human:-1.000000000000000,(Cow:-1.000000000000000,Whale:-1.000000000000000):-1.000000000000000):-1.000000000000000,Seal:-1.000000000000000):-1.000000000000000,(Mouse:-1.000000000000000,Rat:-1.000000000000000):-1.000000000000000):-1.000000000000000,(Platypus:-1.000000000000000,Opossum:-1.000000000000000):-1.000000000000000):-1.000000000000000):-1.000000000000000):-1.000000000000000);";
        // cout << "\nhardcoded tree_string = " << tree_string << endl;
        
    	pllNewickTree *pll_tree = pllNewickParseString(tree_string.c_str());
		assert(pll_tree != NULL);
		pllTreeInitTopologyNewick(ptree->pllInst, pll_tree, PLL_FALSE);
        pllTreeToNewick(ptree->pllInst->tree_string, ptree->pllInst, ptree->pllPartitions, 
                ptree->pllInst->start->back, PLL_TRUE,
			    PLL_TRUE, 0, 0, 0, PLL_SUMMARIZE_LH, 0, 0);  
        string tree_str = string(ptree->pllInst->tree_string);        
        cout << "\nBEFORE eval, string(pllInst->tree_string) = " << tree_str << endl;

        pllNewickParseDestroy(&pll_tree);
		
		ptree->pllInst->bestParsimony = UINT_MAX; // Important because of early termination in evaluateSankoffParsimonyIterativeFastSIMD
		unsigned int pll_score = pllEvaluateParsimonyFast(ptree->pllInst, ptree->pllPartitions, ptree->pllInst->start, PLL_TRUE);
		cout << "\nParsimony score (by PLL kernel) is: " << pll_score << endl;
        
        pllTreeToNewick(ptree->pllInst->tree_string, ptree->pllInst, ptree->pllPartitions, 
                ptree->pllInst->start->back, PLL_TRUE,
			    PLL_TRUE, 0, 0, 0, PLL_SUMMARIZE_LH, 0, 0);             
        tree_str = string(ptree->pllInst->tree_string);        
        cout << "string(pllInst->tree_string) = " << tree_str << endl;

		pllFreeParsimonyDataStructures(ptree->pllInst, ptree->pllPartitions);
        delete ptree;
    }	
}


// Run unit test for parsimony SPR performance:
// ./iqtree2-mpi -s example.phy -t example.phy.treefile --test-pars-spr 1 -redo
// --test-pars-spr <1|2> where 1 is for PLL core, 2 is for IQTREE2 core
void ParsimonyTest::runUnitTestParsimonySPR(Params &params) {
    if(params.unit_test_pars_spr == 2){
        cout << "Test using iqtree2/dev parsimony computation:" << endl;
        Alignment alignment(params.aln_file, params.sequence_type, params.intype, params.model_name);

        // Diep: 2021-05-04
        // The following line is NEW and EXTREMELY IMPORTANT for iqtree2 parsimony.
        // I'm wondering ... Does everyone know?
        alignment.orderPatternByNumChars(PAT_VARIANT);

        IQTree * ptree;
        ptree = new IQTree(&alignment);
        ptree->readTree(params.user_file.c_str(), params.is_rooted); // Read user tree
        ptree->setAlignment(&alignment); // IMPORTANT: Always call setAlignment() after readTree()
        ptree->setParams(&params);
        // ptree->drawTree(cout);

    //    if (ptree->rooted)
    //        ptree->convertToUnrooted();
    //    ptree->initCostMatrix(CM_UNIFORM);
    //    ptree->setParsimonyKernel(params.SSE);

        ptree->initializeAllPartialPars();
        cout << "Parsimony score (by iqtree2 kernel) is: "
                << ptree->computeParsimony()
                << endl;
        delete ptree;
    }else if(params.unit_test_pars_spr == 1){
        cout << "Test using PLL parsimony computation:" << endl;
        Alignment alignment(params.aln_file, params.sequence_type, params.intype, params.model_name);
        // alignment.orderPatternByNumChars(PAT_VARIANT);
        IQTree * ptree;
        ptree = new IQTree(&alignment);
        ptree->readTree(params.user_file.c_str(), params.is_rooted); // Read user tree
        ptree->setAlignment(&alignment); // IMPORTANT: Always call setAlignment() after readTree()
        ptree->setParams(&params);
        
        ptree->initializePLL(params); // Diep: This one modifies topo?
        allocateParsimonyDataStructures(ptree->pllInst, ptree->pllPartitions);

    	// string tree_string = ptree->getTreeString();
        stringstream tree_stream;
        ptree->printTree(tree_stream); // ptree->printTree(tree_stream, WT_SORT_TAXA);
        string tree_string = tree_stream.str();
        cout << "\niqtree2 getTreeString = " << tree_string << endl;
    
        // tree_string = "(LngfishAu:-1.000000000000000,(LngfishSA:-1.000000000000000,LngfishAf:-1.000000000000000):-1.000000000000000,(Frog:-1.000000000000000,((Turtle:-1.000000000000000,((Sphenodon:-1.000000000000000,Lizard:-1.000000000000000):-1.000000000000000,(Crocodile:-1.000000000000000,Bird:-1.000000000000000):-1.000000000000000):-1.000000000000000):-1.000000000000000,((((Human:-1.000000000000000,(Cow:-1.000000000000000,Whale:-1.000000000000000):-1.000000000000000):-1.000000000000000,Seal:-1.000000000000000):-1.000000000000000,(Mouse:-1.000000000000000,Rat:-1.000000000000000):-1.000000000000000):-1.000000000000000,(Platypus:-1.000000000000000,Opossum:-1.000000000000000):-1.000000000000000):-1.000000000000000):-1.000000000000000):-1.000000000000000);";
        // cout << "\nhardcoded tree_string = " << tree_string << endl;
        
    	pllNewickTree *pll_tree = pllNewickParseString(tree_string.c_str());
		assert(pll_tree != NULL);
		pllTreeInitTopologyNewick(ptree->pllInst, pll_tree, PLL_FALSE);
        pllTreeToNewick(ptree->pllInst->tree_string, ptree->pllInst, ptree->pllPartitions, 
                ptree->pllInst->start->back, PLL_TRUE,
			    PLL_TRUE, 0, 0, 0, PLL_SUMMARIZE_LH, 0, 0);  
        string tree_str = string(ptree->pllInst->tree_string);        
        cout << "\nBEFORE eval, string(pllInst->tree_string) = " << tree_str << endl;

        pllNewickParseDestroy(&pll_tree);
		
		ptree->pllInst->bestParsimony = UINT_MAX; // Important because of early termination in evaluateSankoffParsimonyIterativeFastSIMD
		unsigned int pll_score = pllEvaluateParsimonyFast(ptree->pllInst, ptree->pllPartitions, ptree->pllInst->start, PLL_TRUE);
		cout << "\nParsimony score (by PLL kernel) is: " << pll_score << endl;
        
        pllTreeToNewick(ptree->pllInst->tree_string, ptree->pllInst, ptree->pllPartitions, 
                ptree->pllInst->start->back, PLL_TRUE,
			    PLL_TRUE, 0, 0, 0, PLL_SUMMARIZE_LH, 0, 0);             
        tree_str = string(ptree->pllInst->tree_string);        
        cout << "string(pllInst->tree_string) = " << tree_str << endl;
        
		pllFreeParsimonyDataStructures(ptree->pllInst, ptree->pllPartitions);
        delete ptree;
    }	
}

