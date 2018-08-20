
/*
 * Author: Pierre Geurts (p.geurts@ulg.ac.be)
 * (c) 2002-2010
 *
 * (Please do not distribute without permission of the author)
 *
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#include "mex.h"

/* On inclut mon code d'arbre */

#include "tree-model.h"
#include "tree-model.c"
#include "f-table.c"
#include "tree-multiregr.c"

CORETABLE_TYPE *core_table_y;

float getobjy_multiregr_learn_matlab(int obj, int att)
{
  return (float)core_table_y[goal_multiregr[att]*nb_obj_in_core_table+object_mapping[obj]];
}

static void fake_answer(mxArray *plhs[])
{
  plhs[0] = mxCreateDoubleMatrix(0, 0, mxREAL);
}

void set_tree_param(mxArray *treeparam) {
  int nbfield,ifield;
  int nmin=1, bestfirst=0, extratrees=0, adjustdefaultk=1, extratreesk=1, maxnbsplits=5, savepred=1, rf=0;
  float varmin=0.0;
  float a_r=1.0;
  
  char *fname;
  
  if (!mxIsStruct(treeparam))
    mexErrMsgTxt("Fourth argument must be a structure\n");
  
  nbfield=mxGetNumberOfFields(treeparam);
  
  for (ifield=0;ifield<nbfield;ifield++) {
    fname=(char *)mxGetFieldNameByNumber(treeparam,ifield);
    if (strcmp(fname,"nmin")==0) {
      nmin=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"varmin")==0) {
      varmin=(float)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"bestfirst")==0) {
      bestfirst=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"maxnbsplits")==0) {
      maxnbsplits=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"extratrees")==0) {
      extratrees=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"extratreesk")==0) {
      extratreesk=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"rf")==0) {
      rf=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"adjustdefaultk")==0) {
      adjustdefaultk=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"savepred")==0) {      
      savepred=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else
      mexPrintf("Field %s non recognized\n",fname);
  }
  
  init_multiregr_trees(nmin,varmin,a_r,savepred);
  
  set_test_classical();
  
  set_best_first(bestfirst,0,maxnbsplits);

  if (adjustdefaultk==1)
    extratreesk=ceil(sqrt(nb_attributes));
  
  if (extratrees==1) {
    find_a_threshold_num=find_a_threshold_at_random_multiregr;
    find_a_threshold_symb=find_the_best_threshold_symb_multiregr;
  } else {
    find_a_threshold_num=find_the_best_threshold_multiregr;
    find_a_threshold_symb=find_the_best_threshold_symb_multiregr;
  }
  
  if (extratrees==1) {
    find_a_split=find_a_split_at_random_et;
    nb_of_random_tests=extratreesk;
    random_split_score_threshold=10.0;
    rf_k=extratreesk;
  } else if (rf==1) {
    find_a_split=find_the_best_split_among_k;
    nb_of_random_tests=extratreesk;
    random_split_score_threshold=10.0;
    rf_k=extratreesk;
  } else {
    find_a_split=find_the_best_split;
    random_split_score_threshold=10.0;
  }
}

int mart=0;

void set_ensemble_param(mxArray *treeparam) {
  int nbfield,ifield, method;
  int nbterms=1;
  int bootstrap=0;
  float martmu=1.0;

  mart=0;

  char *fname;
  if (!mxIsStruct(treeparam))
    mexErrMsgTxt("Fourth argument must be a structure\n");

  /* analyse les parametres passees en arguments */
  nbfield=mxGetNumberOfFields(treeparam);
  
  for (ifield=0;ifield<nbfield;ifield++) {
    fname=(char *)mxGetFieldNameByNumber(treeparam,ifield);
    if (strcmp(fname,"nbterms")==0) {
      nbterms=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"bootstrap")==0) {
      bootstrap=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"rtparam")==0) {
      set_tree_param(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"mart")==0) {      
      mart=(int)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else if (strcmp(fname,"martmu")==0) {      
      martmu=(float)mxGetScalar(mxGetFieldByNumber(treeparam,0,ifield));
    } else
      mexPrintf("Field %s non recognized\n",fname);
  }

  if (bootstrap==1)
    method=1;
  else if (mart==1)
    method=2;
  else
    method=0;

  set_ensemble_method_parameters(method,nbterms,1,0,martmu);
}

/* cree la structure matlab contenant l'arbre */

/* tres laid pour le moment. Il faudrait peut-etre utiliser
 * la meme structure des deux cotes (surtout si on veut passer
 * le test du cote c)
 */

int save_ensemble_ls_pos_2=0;
int pos_prediction_values=0;
float sum_weight=0.0;

void create_struct_from_tree(int t, mxArray *structvector) {
  mxArray *tmp,*tmp2;
  int nodenumber,start,end,i,j, nb_pred;
  int *ptr_int;
  float *ptr_float;

  /* NODENUMBER */

  start=ltrees[t];
  if (t==current_nb_of_ensemble_terms-1) { /* on est au dernier */
    nodenumber=index_nodes-ltrees[t]+1;
    end=index_nodes+1;
  }else {
    nodenumber=ltrees[t+1]-ltrees[t];
    end=ltrees[t+1];
  }

  tmp=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  *((int *)mxGetPr(tmp))=nodenumber;
  mxSetFieldByNumber(structvector,t,0,tmp);
  
  /* TESTATTRIBUTE */

  tmp=mxCreateNumericMatrix(nodenumber,1,mxINT32_CLASS,mxREAL);
  ptr_int=(int *)mxGetPr(tmp);
  for (i=start; i<end;i++,ptr_int++)
    *ptr_int=tested_attribute[i]+1;
  mxSetFieldByNumber(structvector,t,1,tmp);

  /* TESTTHRESHOLD */

  tmp=mxCreateNumericMatrix(nodenumber,1,mxSINGLE_CLASS,mxREAL);
  ptr_float=(float *)mxGetPr(tmp);
  for (i=start; i<end;i++,ptr_float++)
    *ptr_float=threshold[i].f;
  mxSetFieldByNumber(structvector,t,2,tmp);
  
  /* CHILDREN */

  tmp=mxCreateNumericMatrix(nodenumber,2,mxINT32_CLASS,mxREAL);
  ptr_int=(int *)mxGetPr(tmp);
  for (i=start; i<end;i++,ptr_int++)
    *ptr_int=i+left_successor[i]-start+1;
  for (i=start; i<end;i++,ptr_int++)
    *ptr_int=i+right_successor[i]-start+1;
  mxSetFieldByNumber(structvector,t,3,tmp);
  
  /* LS, OBJECTWEIGHTS */

  tmp=mxCreateNumericMatrix(save_ensemble_ls_size[t],1,mxINT32_CLASS,mxREAL);
  tmp2=mxCreateNumericMatrix(save_ensemble_ls_size[t],1,mxSINGLE_CLASS,mxREAL);
  ptr_int=(int *)mxGetPr(tmp);
  ptr_float=(float *)mxGetPr(tmp2);
  for (i=0; i<save_ensemble_ls_size[t];i++,ptr_int++,ptr_float++) {
    *ptr_int=save_ensemble_ls_vector[save_ensemble_ls_pos_2]+1;
    *ptr_float=save_ensemble_ls_weight[save_ensemble_ls_pos_2];
    save_ensemble_ls_pos_2++;
  }
  mxSetFieldByNumber(structvector,t,4,tmp);
  mxSetFieldByNumber(structvector,t,5,tmp2);

  /* WEIGHT */

  tmp=mxCreateNumericMatrix(1,1,mxSINGLE_CLASS,mxREAL);
  if (average_predictions_ltrees)
    *((float *)mxGetPr(tmp))=ltrees_weight[t]/sum_weight;
  else
    *((float *)mxGetPr(tmp))=ltrees_weight[t];
  mxSetFieldByNumber(structvector,t,6,tmp);

  /* SAVEPRED */
  tmp=mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
  *((int *)mxGetPr(tmp))=multiregr_savepred;
  mxSetFieldByNumber(structvector,t,7,tmp);

  if (multiregr_savepred) {
    /* INDEXPREDICTION */
    tmp=mxCreateNumericMatrix(nodenumber,1,mxINT32_CLASS,mxREAL);
    ptr_int=(int *)mxGetPr(tmp);
    for (i=start; i<end;i++,ptr_int++) {
      if (prediction[i]<0)
	*ptr_int=0;
      else
	*ptr_int=prediction[i]-pos_prediction_values+1;
    }
    mxSetFieldByNumber(structvector,t,8,tmp);
    
    /* PREDICTIONS */
    nb_pred=(nodenumber+1)/2;
    tmp=mxCreateNumericMatrix(nb_pred,nb_goal_multiregr,mxSINGLE_CLASS,mxREAL);
    ptr_float=(float *)mxGetPr(tmp);
    for (i=0;i<nb_goal_multiregr;i++) {
      for (j=pos_prediction_values;j<pos_prediction_values+nb_pred;j++) {
	*ptr_float=prediction_values[j][i];
	ptr_float++;
      }
    }
    mxSetFieldByNumber(structvector,t,9,tmp);
    
    pos_prediction_values+=nb_pred;
    
  }

  /* NODESIZE */
  tmp=mxCreateNumericMatrix(nodenumber,1,mxSINGLE_CLASS,mxREAL);
  ptr_float=(float *)mxGetPr(tmp);
  for (i=start; i<end;i++,ptr_float++)
    *ptr_float=node_size[i];
  mxSetFieldByNumber(structvector,t,10,tmp);
  
}

mxArray *create_struct_from_current_ensemble() {
  mxArray *structvector, *structens, *tmp;
  int t, *ptr_int,i;
  char const* fieldnames[]={"nodenumber","testattribute","testthreshold","children","ls","objectweights","weight","savepred","indexprediction","predictions","nodesize"};
  char const* fieldnamesens[]={"ls","trees"};

  save_ensemble_ls_pos_2=0;
  pos_prediction_values=0;
  sum_weight=0.0;
  if (average_predictions_ltrees) {
    for (t=0; t<current_nb_of_ensemble_terms; t++)
      sum_weight+=ltrees_weight[t];
  }

  structvector=mxCreateStructMatrix(current_nb_of_ensemble_terms,1,11,fieldnames);
  
  for (t=0; t<current_nb_of_ensemble_terms; t++) {
    /* cree la structure pour chaque arbre */
    create_struct_from_tree(t,structvector);
  }

  structens=mxCreateStructMatrix(1,1,2,fieldnamesens);

  tmp=mxCreateNumericMatrix(global_learning_set_size,1,mxINT32_CLASS,mxREAL);
  ptr_int=(int *)mxGetPr(tmp);
  for (i=0; i<global_learning_set_size;i++,ptr_int++)
    *ptr_int=object_mapping[i]+1;

  mxSetFieldByNumber(structens,0,0,tmp);
  mxSetFieldByNumber(structens,0,1,structvector);

  return structens;
}

/* test l'ensemble sur un test set */

mxArray *compute_ts_predictions_matlab_multiregr(int nts,float *xts) {
  mxArray *table;
  float *ptr_float, *tmp;
  float *saved_core_table=core_table;
  int saved_nb_obj_in_core_table=nb_obj_in_core_table;
  int i,j;
  void (*get_ens_pred)(int obj, float *vect_pred);

  /* savepred */
  if (multiregr_savepred==0) {
    /* calcul les stats aux feuilles */
    compute_node_subset_current_ensemble_weighted(NULL);
    get_ens_pred=get_multiregr_ens_pred;
  } else {
    get_ens_pred=get_multiregr_ens_pred_savepred;
  }
  /* on change la core table */
  core_table=xts;
  nb_obj_in_core_table=nts;

  /* on cree la table */
  table=mxCreateNumericMatrix(nts,nb_goal_multiregr,mxSINGLE_CLASS,mxREAL);
  ptr_float=(float *)mxGetPr(table);

  /* on cree un vecteur pour contenir la prediction */
  tmp=(float *)mxMalloc((size_t)nb_goal_multiregr*sizeof(float));

  /* on teste */
  for (i=0; i<nts; i++) {
    get_ens_pred(i,tmp);
    for (j=0; j<nb_goal_multiregr; j++) {
      ptr_float[j*nts+i]=tmp[j];
    }
  }
  
  mxFree(tmp);

  /* on remet la coretable (pour les variables importances) */
  core_table=saved_core_table;
  nb_obj_in_core_table=saved_nb_obj_in_core_table;
  
  return table;
}

mxArray *compute_variable_imp_matlab() {
  mxArray *varimp;
  double *varimpptr;
  int  i;

  /* on retrie les variables */
  for (i=0; i<nb_attributes;i++)
    attribute_vector[i]=i;

  varimp=mxCreateNumericMatrix(nb_attributes,nb_goal_multiregr,mxDOUBLE_CLASS,mxREAL);
  varimpptr=(double *)mxGetPr(varimp);
  for (i=0;i<nb_attributes*nb_goal_multiregr; i++){
    varimpptr[i]=0.0;
  }
  compute_ltrees_variable_importance_multiregr_separate(varimpptr,-1);    
  
  return varimp;
}

/* Interface function of matlab */

/* now assume prhs[0]: label prhs[1]: features */

void mexFunction(int nlhs, mxArray *plhs[],
		 int nrhs, const mxArray *prhs[])
{
  mxArray *xData,*yData, *wData, *lsData, *treeparamData, *goaltypematlab, *attdes;
  float *xValues, *yValues;
  int i,j,pos,verbose,length_ls_vector,nbweights,maxnbnodes, *ptr_int;
  
  if ((nrhs<6)||(nrhs>7)) {
    mexErrMsgTxt("\n\nToo many or too few arguments\nThis function takes the following arguments\n  - X: a single float table with input data\n  - Y: a single float table with output data\n  - ls: the index of the objects in X and Y that are using for growing the tree\n  - w: a vector of object weights\n  - param: a structure that contains the parameter of the method\n  (- XTS: a single float table with input data for the test case)\n  - verbose: 0 to disable outputs during training\n\nThe function returns:\n  (- YTS: output predictions for the test examples if XTS is provided)\n  - a matlab structure containing the tree\n  - Variable importance as measured by the function (not yet implemented)\n\nBe carefull, argument types are not checked\n ");
  }

  xData=(mxArray *)prhs[0];
  yData=(mxArray *)prhs[1];
  lsData=(mxArray *)prhs[2];
  wData=(mxArray *)prhs[3];
  treeparamData=(mxArray *)prhs[4];
  if (nrhs==6)
    verbose=(int)mxGetScalar(prhs[5]);
  else
    verbose=(int)mxGetScalar(prhs[6]);

  set_print_result(verbose,0);

  /* initialise les entrees */
  nb_attributes=mxGetN(xData);
  nb_obj_in_core_table=mxGetM(xData);  
  core_table=(CORETABLE_TYPE *)mxGetPr(xData);

  attribute_descriptors=(int *)mxMalloc((size_t)nb_attributes*sizeof(int));

  /* */

  length_attribute_descriptors=nb_attributes;
  
  attdes=mxGetField(treeparamData,0,"attdesc");
  if (attdes==NULL) {
    for (i=0;i<nb_attributes;i++)
      attribute_descriptors[i]=0; /* on suppose tous les att numeriques, on pourrait facilement traiter des symboliques */  
  } else {
    ptr_int=(int *)mxGetPr(attdes);
    for (i=0; i<nb_attributes; i++,ptr_int++)
      attribute_descriptors[i]=(*ptr_int);
  } 

  nb_classes=0;
  attribute_vector=mxMalloc((size_t)nb_attributes*sizeof(int));
  for (i=0; i<nb_attributes;i++)
    attribute_vector[i]=i;

  /* on construit le LS et le vecteur de poids (init_learning_set) */

  nb_obj_in_core_table=mxGetM(xData);
  length_ls_vector=mxGetN(lsData);
  global_learning_set_size=length_ls_vector;
  current_learning_set_size=length_ls_vector;
  object_mapping=(int *)mxMalloc((size_t)length_ls_vector*sizeof(int));
  ptr_int=(int *)mxGetPr(lsData);
  for (i=0; i<length_ls_vector; i++,ptr_int++)
    object_mapping[i]=(*ptr_int)-1;

  current_learning_set=(int *)mxMalloc((size_t)length_ls_vector*sizeof(int));
  for (i=0; i<length_ls_vector; i++)
    current_learning_set[i]=i;

  nbweights=mxGetN(wData);

  if (nbweights==0) {
    object_weight=(SCORE_TYPE *)mxMalloc((size_t)length_ls_vector*sizeof(SCORE_TYPE));
    for (i=0;i<length_ls_vector;i++)
      object_weight[i]=1.0;
  } else if (nbweights==length_ls_vector) {
    object_weight=(double *)mxGetPr(wData);
  } else
    mexErrMsgTxt("wrong size of the weight vector");


  /* initialise un probleme de regression multiple */
  nb_goal_multiregr=mxGetN(yData);
  core_table_y=(CORETABLE_TYPE *)mxGetPr(yData);
  goal_multiregr=mxMalloc((size_t)nb_goal_multiregr*sizeof(int));
  for (i=0;i<nb_goal_multiregr; i++)
    goal_multiregr[i]=i;
  getobjy_multiregr_learn=getobjy_multiregr_learn_matlab;
  
  /* on fixe les parametres des arbres */
  set_ensemble_param(treeparamData);

  init_save_ensemble_ls(1);

  /* alloue la memoire pour les arbres */
  maxnbnodes=number_of_ensemble_terms*(best_first*(2*best_first_max_nb_tests+1)+(1-best_first)*(2*length_ls_vector-1));
  allocate_tree_tables(maxnbnodes,ceil((maxnbnodes+number_of_ensemble_terms)/2),multiregr_savepred*nb_goal_multiregr,0);
  allocate_multiregr_table_score(nb_goal_multiregr);
    
  clean_all_trees();

  /* construction de l'ensemble d'arbres */
  
  build_one_tree_ensemble(NULL, 0);

  
  /* on retourne le tout */

  if (nlhs>0) {
    if (nrhs==7) {
      if (verbose==1) {
	mexPrintf("Test...");
      }

      plhs[0]=compute_ts_predictions_matlab_multiregr(mxGetM(prhs[5]),(float *)mxGetPr(prhs[5]));

      if (verbose==1) {
	mexPrintf("\n");
      }
      if (nlhs>1) {
	if (verbose==1) {
	  mexPrintf("Calcul de l'importance des variables...");
	}
	plhs[1]=compute_variable_imp_matlab();
	if (verbose==1) {
	  mexPrintf("\n");
	}
      }
      if (nlhs>2) {/* on renvoit aussi la structure */
	if (verbose==1) {
	  mexPrintf("Creation de la structure...");
	}
	plhs[2]=create_struct_from_current_ensemble();
	if (verbose==1) {
	  mexPrintf("\n");
	}
      }      
    } else {
      if (verbose==1) {
	mexPrintf("Creation de la structure...");
      }
      plhs[0]=create_struct_from_current_ensemble();
      if (verbose==1) {
	mexPrintf("\n");
      }
      if (nlhs>1) {
	if (verbose==1) {
	  mexPrintf("Calcul de l'importance des variables...");
	}
	plhs[1]=compute_variable_imp_matlab();
	if (verbose==1) {
	  mexPrintf("\n");
	}
      }
    }
  }

  /* on libere la memoire */
  if (verbose==1)
    mexPrintf("Liberation de la memoire\n");
  mxFree(attribute_descriptors);
  mxFree(attribute_vector);
  mxFree(goal_multiregr);
  mxFree(object_mapping);
  mxFree(current_learning_set);

  if (nbweights==0)
    mxFree(object_weight);
  free_tree_tables();
  MyFree((int *)save_ensemble_ls_vector);
  save_ensemble_ls_vector=NULL;
  MyFree((float *)save_ensemble_ls_weight);
  save_ensemble_ls_weight=NULL;

  return;
}
