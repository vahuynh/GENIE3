
/*
 * Author: Pierre Geurts (p.geurts@ulg.ac.be)
 * (c) 2002-2010
 *
 * (Please do not distribute without permission of the author)
 *
 */


/******************************/
/* ALGORITHME DE CONSTRUCTION */
/******************************/

int multiregr_savepred=0;

float (*getobjy_multiregr_learn)(int obj, int i); /* permet de separer la table des entrees de la table des sorties */

int allocate_multiregr_table_score(int nb_goal) {

  table_score_multiregr=allocate_table_score_type(3,2*nb_goal+1);
  if (table_score_multiregr==NULL) {
    free_tree_tables();
    return 0;
  }
  
  table_score_symb_multiregr=allocate_table_score_type(MAX_NUMBER_OF_SYMBOLIC_VALUES,2*nb_goal+1);
  if (table_score_symb_multiregr==NULL) {
    free_tree_tables();
    return 0;
  }
  return 1;
}

void free_multiregr_table_score() {
  int i;

  if (table_score_multiregr!=NULL) {
    for (i=0;i<3;i++) {
      if (table_score_multiregr[i]!=NULL)
	MyFree((SCORE_TYPE *)table_score_multiregr[i]);
    }
    MyFree((SCORE_TYPE **)table_score_multiregr);
    table_score_multiregr=NULL;
  }
  if (table_score_symb_multiregr!=NULL) {
    for (i=0;i<MAX_NUMBER_OF_SYMBOLIC_VALUES;i++) {
      if (table_score_symb_multiregr[i]!=NULL)
	MyFree((SCORE_TYPE *)table_score_symb_multiregr[i]);
    }
    MyFree((SCORE_TYPE **)table_score_symb_multiregr);
    table_score_symb_multiregr=NULL;
  }
}

int *goal_multiregr;
int nb_goal_multiregr;

float getobjy_multiregr_learn_true(int obj, int i) {
  return getattval(object_mapping[obj],goal_multiregr[i]);
}

/* calcule les stats sur les sorties */

void summarize_vector_multiregr(int *vector, int start, int end) {
  int i,j;
  SCORE_TYPE w;
  
  
  for (i=0; i<2*nb_goal_multiregr+1; i++) {
    table_score_multiregr[0][i]=0.0;
  }

  for (i=start; i<= end; i++) {
    w=object_weight[vector[i]];
    table_score_multiregr[0][0]+=w;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(vector[i],j);
      table_score_multiregr[0][2*j+1]+=w*y;
      table_score_multiregr[0][2*j+2]+=w*y*y;
    }
  }

  /*
    printf(" N=%e\n",table_score_multiregr[0][0]);
    for (j=0; j<nb_goal_multiregr; j++) {
    printf("%d %e %e\n", j, table_score_multiregr[0][2*j+1], table_score_multiregr[0][2*j+2]);
    }
    fflush(stdout);
  */
}

/* fonction de repropagation des objets dans l'arbre */

int node_subset_desc_size=0;
typedef int nsd_type[2];

nsd_type *node_subset_desc=NULL; /* [0] is the number of cases in that node
				  * [1] is the index in the node_object_list
				  */

void free_node_subset_desc() {
  if (node_subset_desc_size>0)
    MyFree((nsd_type *)node_subset_desc);
}

void allocate_node_subset_desc(int size) {
  if (node_subset_desc_size<size) {
    free_node_subset_desc();
    node_subset_desc=(nsd_type *)MyMalloc((size_t)size*2*sizeof(int));
    if (!node_subset_desc) {
      printf("impossible d'allouer de la memoire dans allocate_node_subset_desc\n");
      exit(0);
    }
    node_subset_desc_size=size;
  }
}

int *node_object_list;
int node_object_list_size=0;

/* pour renvoyer les resultats au lisp */
float *object_weight2;

float *node_object_weight;
int node_object_allocatedp=0;


/* SCORE_TYPE v_tot; */

/* le score est simplement la somme des réductions de variance */

SCORE_TYPE compute_multiregr_score_from_table() {
  int i;
  SCORE_TYPE y_tot_var, n_tot_var;

  y_tot_var=0.0;
  n_tot_var=0.0;

  table_score_multiregr[2][0]=table_score_multiregr[0][0]-table_score_multiregr[1][0];
  
  for (i=0; i<nb_goal_multiregr; i++) {
    table_score_multiregr[2][2*i+1]=table_score_multiregr[0][2*i+1]-table_score_multiregr[1][2*i+1];
    table_score_multiregr[2][2*i+2]=table_score_multiregr[0][2*i+2]-table_score_multiregr[1][2*i+2];
    
    y_tot_var+=fabs(table_score_multiregr[1][2*i+2]-(table_score_multiregr[1][2*i+1]*table_score_multiregr[1][2*i+1])/table_score_multiregr[1][0]);
    n_tot_var+=fabs(table_score_multiregr[2][2*i+2]-(table_score_multiregr[2][2*i+1]*table_score_multiregr[2][2*i+1])/table_score_multiregr[2][0]);
  }
  
  info=v_tot-(y_tot_var+n_tot_var);

  return (info/v_tot); /* le score est compris entre 0 et 1. La
			  normalisation n'est pas nécessaire */
}

SCORE_TYPE compute_multiregr_score_from_table_min() {
  int i;
  SCORE_TYPE y_tot_var, n_tot_var, v_tot_cur, y_tot_var_cur, n_tot_var_cur, min_score, min_info, score_cur;

  /* y_tot_var=0.0; */
  /* n_tot_var=0.0; */
  v_tot_cur=0.0;
  y_tot_var_cur=0.0;
  n_tot_var_cur=0.0;
  min_score=10.0;
  min_info=0.0;

  table_score_multiregr[2][0]=table_score_multiregr[0][0]-table_score_multiregr[1][0];
  
  for (i=0; i<nb_goal_multiregr; i++) {
    table_score_multiregr[2][2*i+1]=table_score_multiregr[0][2*i+1]-table_score_multiregr[1][2*i+1];
    table_score_multiregr[2][2*i+2]=table_score_multiregr[0][2*i+2]-table_score_multiregr[1][2*i+2];
    
    v_tot_cur=table_score_multiregr[0][2*i+2]-(table_score_multiregr[0][2*i+1]*table_score_multiregr[0][2*i+1])/table_score_multiregr[0][0];
    
    y_tot_var_cur=fabs(table_score_multiregr[1][2*i+2]-(table_score_multiregr[1][2*i+1]*table_score_multiregr[1][2*i+1])/table_score_multiregr[1][0]);
    n_tot_var_cur=fabs(table_score_multiregr[2][2*i+2]-(table_score_multiregr[2][2*i+1]*table_score_multiregr[2][2*i+1])/table_score_multiregr[2][0]);
    score_cur=(v_tot_cur-(y_tot_var_cur+n_tot_var_cur))/v_tot_cur;

    if (score_cur<=min_score) {
      min_score=score_cur;
      min_info=(v_tot_cur-(y_tot_var_cur+n_tot_var));
    }
  }
  
  info=min_info;

  return min_score; /* le score est compris entre 0 et 1. La
			  normalisation n'est pas nécessaire */
}

/* stop splitting criterio pour la distance */

/* SCORE_TYPE v_min=0.0; */
/* verifier ce calcul) */

int stop_splitting_criterio_multiregr() {
  int i;

  v_tot=0.0;
  for (i=0; i<nb_goal_multiregr; i++) {
    /* printf("var %d %e\n", i, table_score_multiregr[0][2*i+2]-(table_score_multiregr[0][2*i+1]*table_score_multiregr[0][2*i+1])/table_score_multiregr[0][0]);*/
    v_tot+=table_score_multiregr[0][2*i+2]-(table_score_multiregr[0][2*i+1]*table_score_multiregr[0][2*i+1])/table_score_multiregr[0][0];
  }

  if (print_detailled_result) {
    printf("N=%e, v_tot=%e\n",table_score_multiregr[0][0], v_tot);
    fflush(stdout);
  }

  return ((v_tot/table_score_multiregr[0][0])<=v_min);
}

int not_significant_test_multiregr() {
  PRECISION f_stat;

  /* sans doute tout à fait incorrect */
  if (best_threshold_score>=0.0) {
    f_stat=(1.0/(1.0-best_threshold_score))*(table_score_multiregr[0][0]-2)/(table_score_multiregr[0][0]-1);
    
    if (print_detailled_result) {
      printf("N=%e, Total variance=%e, Variance reduction=%e, score=%e\n", table_score_multiregr[0][0], v_tot, best_threshold_score*v_tot, best_threshold_score);
      printf("F-stat=%e, F-thres=%e\n", f_stat, f_test_threshold(1,table_score_multiregr[0][0]-2, alfa_rt));
      fflush(stdout);
    }
    return (f_stat<f_test_threshold(1,table_score_multiregr[0][0]-2, alfa_rt));
  } else {
    return 1;
  }
}

SCORE_TYPE best_first_node_priority_function_multiregr() {
  int i;
  SCORE_TYPE pr=0.0;

  for (i=0; i<nb_goal_multiregr; i++) {
    pr+=table_score_multiregr[0][2*i+2]-(table_score_multiregr[0][2*i+1]*table_score_multiregr[0][2*i+1])/table_score_multiregr[0][0];
  }
  
  return pr;
}

/* on n'enregistre rien pour le moment, la prédiction sera calculée au moment du test. Ca evite de retenir un grand vecteur 
 * pour chaque feuille (mais ça augmente les temps de calcul pour le test)
 */

int make_leaf_prediction_multiregr_nosave() {
  return -1;
}

int make_leaf_prediction_multiregr_savepred() {
  int i;
  index_prediction_values++;

  for (i=0;i<nb_goal_multiregr; i++) {
    prediction_values[index_prediction_values][i]=(float)(table_score_multiregr[0][2*i+1]/table_score_multiregr[0][0]);
  }
  return index_prediction_values;
}

/* recherche le meilleur seuil pour la distance */

/* attribut numérique */

void find_the_best_threshold_multiregr(int att, int *ls_vector, int start, int end) {
  float old_val, new_val;
  SCORE_TYPE best_score=-1.0, best_info, current_score, w;
  float best_threshold;
  int st=start,i;

  /* initialisation de la table */
  table_score_multiregr[1][0]=0.0;
  for (i=0; i<nb_goal_multiregr; i++) {
    table_score_multiregr[1][2*i+1]=0.0;
    table_score_multiregr[1][2*i+2]=0.0;
  }
  
  /* on trie l'ensemble selon l'attribut */
  sort_ls_vector(ls_vector, start, end, att);

  /* on parcourt toutes les valeurs de seuils possibles */
  old_val=getattval(object_mapping[ls_vector[start]],att);
  
  for(st=start; st<end; st++) {
    w=object_weight[ls_vector[st]];
    table_score_multiregr[1][0]+=w;
    for (i=0; i<nb_goal_multiregr; i++) {
      float y=getobjy_multiregr_learn(ls_vector[st],i);
      table_score_multiregr[1][2*i+1]+=w*y;
      table_score_multiregr[1][2*i+2]+=w*(y*y);
    }

    if ((new_val=getattval(object_mapping[ls_vector[st+1]],att))!=old_val) { /* un nouveau seuil a considerer */
      
      current_score=compute_score_from_table();
      if (current_score>best_score) {
	best_score=current_score;
	best_info=info;
	best_threshold=(old_val+new_val)/2.0;
	if (old_val>=best_threshold) /* problem d'arrondi */
	  best_threshold=new_val;
      }
      old_val=new_val;
    }
  }
  if (best_score>=0.0) {
    current_threshold.f=best_threshold;
    current_threshold_score=best_score;
    current_threshold_info=best_info;
  } else {
    current_threshold_score=-1.0;
  }
}

/* attribut symbolique */

/* attribut symbolique */

void summarize_symb_att_multiregr(int att, int *vector, int start, int end) {
  int i,j;

  /* set to zero */
  for (i=0; i<NB_VAL(att); i++) {
    table_score_symb_multiregr[i][0]=0.0;
    for (j=0; j<nb_goal_multiregr; j++) {
      table_score_symb_multiregr[i][2*j+1]=0.0;
      table_score_symb_multiregr[i][2*j+2]=0.0;
    }
  }
  
  /* fill the table with frequency */
  for (i=start; i<=end; i++) {
    SCORE_TYPE w=object_weight[vector[i]];
    int v=(int)getattval(object_mapping[vector[i]],att);
    table_score_symb_multiregr[v][0]+=w;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(vector[i],j);
      table_score_symb_multiregr[v][2*j+1]+=w*y;
      table_score_symb_multiregr[v][2*j+2]+=w*(y*y);
    }
  }
}

/* ATTENTION version super inefficace et naive */

void find_the_best_threshold_symb_multiregr(int att, int *ls_vector, int start, int end) {
  int i, v, eff_v;
  int nb_val=NB_VAL(att);
  int nb_val_ls=0;
  union threshold_type current_subset, best_subset;
  SCORE_TYPE best_score=-1.0, best_info, current_score;

  /* on precalcule la table avec tous les frequences pour toutes les classes */
  summarize_symb_att_multiregr(att, ls_vector, start, end);

  /* Check that all elements do not belong to the same class */
  for (i=0; i<nb_val; i++) {
    if (table_score_symb_multiregr[i][0]!=0) 
      nb_val_ls++;
  }
  if (nb_val_ls==1) { /* all objects have the same value of this attribute */
    current_threshold_score=-1.0;
    return;
  }

  init_threshold_type(&current_subset);
  add1_threshold_type(&current_subset);

  do {
    /* fill the table score according to the current subset */
    table_score_multiregr[1][0]=0.0;
    for (i=0; i<nb_goal_multiregr; i++) {
      table_score_multiregr[1][2*i+1]=0.0;
      table_score_multiregr[1][2*i+2]=0.0;
    }
      
    eff_v=0;
    for (v=0; v<nb_val; v++) {
      if (table_score_symb_multiregr[v][0]!=0.0) {
	/* check bit eff_v in current_subset */
	if (BITN(current_subset, eff_v)) {
	  table_score_multiregr[1][0]+=table_score_symb_multiregr[v][0];
	  for (i=0; i<nb_goal_multiregr; i++) {
	    table_score_multiregr[1][2*i+1]+=table_score_symb_multiregr[v][2*i+1];
	    table_score_multiregr[1][2*i+2]+=table_score_symb_multiregr[v][2*i+2];
	  }
	}
	eff_v++;
      }
    }
 
    /* compute the score */
    current_score=compute_score_from_table();

    if (current_score>best_score) {
      best_score=current_score;
      best_info=info;
      best_subset=current_subset;
    }
    add1_threshold_type(&current_subset);
  } while (!(BITN(current_subset,(nb_val_ls-1))));
	 
  if (best_score>=0.0) {
    current_threshold_score=best_score;
    current_threshold_info=best_info;
    /* translate current_subset into a proper subset */
    init_threshold_type(&current_threshold);
    eff_v=0;
    for (v=0; v<nb_val; v++) {
      if (table_score_symb_multiregr[v][0]!=0.0) {
	if (BITN(best_subset, eff_v)) {
	  SET_BITN(current_threshold,v);
	}
	eff_v++;
      }
    }    
  } else {
    current_threshold_score=-1.0;
  }
}

/* version extra-trees */

/* attribut numérique */

void find_a_threshold_at_random_multiregr(int att, int *ls_vector, int start, int end) {
  int i, j;
  float min=getattval(object_mapping[ls_vector[start]],att);
  float max=min;
  SCORE_TYPE w;
  
  current_threshold_score=-1.0;

  /* calcule les stats sur l'attribut */
  for (i=start+1; i<=end; i++) {
    float val=getattval(object_mapping[ls_vector[i]],att);
    if (val<min)
      min=val;
    else if (val>max)
      max=val;
  }
  
  if (min==max) { /* toutes les valeurs sont egales */
    return;
  }
  
  /* tirage du seuil (uniformément entre min et max) */
  current_threshold.f=max-(max-min)*get_random_float();
  
  /* calcul du score */
  table_score_multiregr[1][0]=0.0;
  for (i=0; i<nb_goal_multiregr; i++) {
    table_score_multiregr[1][2*i+1]=0.0;
    table_score_multiregr[1][2*i+2]=0.0;
  }

  for (i=start; i<=end; i++) {
    if (getattval(object_mapping[ls_vector[i]],att)<current_threshold.f) {
      w=object_weight[ls_vector[i]];
      table_score_multiregr[1][0]+=w;
      for (j=0; j<nb_goal_multiregr; j++) {
	float y=getobjy_multiregr_learn(ls_vector[i],j);
	table_score_multiregr[1][2*j+1]+=w*y;
	table_score_multiregr[1][2*j+2]+=w*y*y;
      }
    }
  }
  current_threshold_score=compute_score_from_table();
}

void find_a_threshold_symb_at_random_multiregr(int att, int *ls_vector, int start, int end) {
  int i, v, eff_v;
  int nb_val=NB_VAL(att);
  int nb_val_ls=0;
  union threshold_type random_subset;

  /* on precalcule la table avec tous les frequences pour toutes les classes */
  summarize_symb_att_multiregr(att, ls_vector, start, end);
  
  /* Check that all elements do not belong to the same class */
  for (i=0; i<nb_val; i++) {
    if (table_score_symb_multiregr[i][0]!=0) 
      nb_val_ls++;
  }
  if (nb_val_ls==1) { /* all objects have the same value of this attribute */
    current_threshold_score=-1.0;
    return;
  } 

  /* choisit un subset aleatoirement (entre 1 et 2**nb_val_ls-2 inclus) */
  get_random_subset_non_empty_full(&random_subset, nb_val_ls);

  /* on l'insere dans un subset aleatoire sur toutes les valeurs */
  init_threshold_type(&current_threshold);
  eff_v=0;
  for (v=0; v<nb_val; v++) {
    if (table_score_symb_multiregr[v][0]!=0.0) {
      if (BITN(random_subset, eff_v)) {
	SET_BITN(current_threshold,v);
      }
      eff_v++;
    } else {
      if (get_random_float()<0.5) { /* on flippe aleatoirement le bit */
	SET_BITN(current_threshold,v);
      }
    }
  }

  /* compute the score */
  /* fill the table score according to the current subset */
  table_score_multiregr[1][0]=0.0;
  for (i=0; i<nb_goal_multiregr; i++) {
    table_score_multiregr[1][2*i+1]=0.0;
    table_score_multiregr[1][2*i+2]=0.0;
  }
    
  eff_v=0;
  for (v=0; v<nb_val; v++) {
    if (table_score_symb_multiregr[v][0]!=0.0) {
      /* check bit eff_v in current_subset */
      if (BITN(random_subset, eff_v)) {
	table_score_multiregr[1][0]+=table_score_symb_multiregr[v][0];
	for (i=0; i<nb_goal_multiregr; i++) {
	  table_score_multiregr[1][2*i+1]+=table_score_symb_multiregr[v][2*i+1];
	  table_score_multiregr[1][2*i+2]+=table_score_symb_multiregr[v][2*i+2];
	}
      }
      eff_v++;
    }
  }

  /* compute the score */
  current_threshold_score=compute_score_from_table();
}

/**************************/
/* INTERFACE AVEC LE LISP */
/**************************/

void init_multiregr_trees(int n_min, float vmin, float ar, int savepred) {
  min_node_size=n_min;
  v_min=vmin;
  alfa_rt=ar;
  summarize_vector=summarize_vector_multiregr;
  
  multiregr_savepred=savepred;

  if (savepred==1)
    make_leaf_prediction=make_leaf_prediction_multiregr_savepred;
  else
    make_leaf_prediction=make_leaf_prediction_multiregr_nosave;

  stop_splitting_criterio=stop_splitting_criterio_multiregr;
  best_first_node_priority_function=best_first_node_priority_function_multiregr;
  not_significant_test=not_significant_test_multiregr;
  compute_score_from_table=compute_multiregr_score_from_table;
}

/****************************/
/* IMPORTANCE DES VARIABLES */
/****************************/

/* code separé par rapport aux arbres classiques mais les deux pourraient être combinés sans problème */
/* pour le moment, ts_vector contient des index dans object_mapping (de 0 à N)
 */

/* il faudrait l'allouer de manière dynamique */
int attribute_position[MAX_NUMBER_OF_ATTRIBUTES];

SCORE_TYPE get_vi_multiregr(int *ts_vector, int start, int end, int borne) {
  int i,j, st;
  SCORE_TYPE current_score;

  /* summarize_vector */
  for (i=0; i<2*nb_goal_multiregr+1; i++) {
    table_score_multiregr[0][i]=0.0;
    table_score_multiregr[1][i]=0.0;
  }

  for (i=start; i<=end; i++) {
    table_score_multiregr[0][0]++;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(ts_vector[i],j);
      table_score_multiregr[0][2*j+1]+=y;
      table_score_multiregr[0][2*j+2]+=y*y;
    }
  }
  
  /* calcule v_tot en fait */
  stop_splitting_criterio_multiregr();
  
  if ((start>=borne)||(borne>end)) {
    return 0.0;
  }

  /* fill the table_score, il y a sûrement une manière plus efficace que de le faire
   * objet par objet */

  for (i=start; i<borne; i++) {
    table_score_multiregr[1][0]++;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(ts_vector[i],j);
      table_score_multiregr[1][2*j+1]+=y;
      table_score_multiregr[1][2*j+2]+=(y*y);
    }
  }

  /* compute the score */
  current_score=compute_multiregr_score_from_table();
  
  return (info);
}


SCORE_TYPE compute_one_tree_variable_importance_multiregr(int tree, int *ts_vector, int length_ts_vector, float weight,
							  SCORE_TYPE *attribute_importance, int obj) {
  SCORE_TYPE total_var=-1.0;
  
  index_stack_open_nodes++;
  stack_open_nodes[index_stack_open_nodes][0]=tree;
  stack_open_nodes[index_stack_open_nodes][1]=0;
  stack_open_nodes[index_stack_open_nodes][2]=length_ts_vector-1;
  
  while (index_stack_open_nodes>=0) {
    int node=stack_open_nodes[index_stack_open_nodes][0];
    int start=stack_open_nodes[index_stack_open_nodes][1];
    int end=stack_open_nodes[index_stack_open_nodes][2];
    int node_size=end-start+1;
    float vi;

    if ((left_successor[node]==-1)||(node_size==1)) {
      index_stack_open_nodes--;
      /*
	printf("term start=%d, end=%d, tested attribute=%d isop=%d\n", start, end, tested_attribute[node], index_stack_open_nodes);
	fflush(stdout);
      */
    } else {
      /* separation */
      int borne=separate_ls_vector_local(tested_attribute[node], threshold[node], ts_vector, start, end);

      /*
	printf("current_node=%d, left=%d, right=%d\n",node, left_successor[node], right_successor[node]);
	printf("start=%d, borne=%d, end=%d, tested attribute=%d isop=%d\n", start, borne, end, tested_attribute[node], index_stack_open_nodes);
	fflush(stdout);
      */

      /* calcul de l'importance (seulement si borne OK */
      vi=get_vi_multiregr(ts_vector, start, end, borne);
      
      /*
      printf("start=%d borne=%d end=%d vi=%f\n",start,borne,end, vi);
      fflush(stdout);
      */
      
      if (node==tree) {
	  total_var=v_tot;
      }
	/* mis a jour du vecteur */
      attribute_importance[attribute_position[tested_attribute[node]]]+=(weight*vi);

      /* left and right successors are put on the stack */
      index_stack_open_nodes--;
      if (obj<0) {
	if (start<borne) {
	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=node+left_successor[node];
	  stack_open_nodes[index_stack_open_nodes][1]=start;
	  stack_open_nodes[index_stack_open_nodes][2]=borne-1;
	}
	if (borne<=end) {
	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=node+right_successor[node];
	  stack_open_nodes[index_stack_open_nodes][1]=borne;
	  stack_open_nodes[index_stack_open_nodes][2]=end;
	}
      } else {
	if (check_test(tested_attribute[node],getattval(obj,tested_attribute[node]),threshold[node])) {
	  /* on met le gauche */
	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=node+left_successor[node];
	  stack_open_nodes[index_stack_open_nodes][1]=start;
	  stack_open_nodes[index_stack_open_nodes][2]=borne-1;
	} else {
	  /* on met le droit */
	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=node+right_successor[node];
	  stack_open_nodes[index_stack_open_nodes][1]=borne;
	  stack_open_nodes[index_stack_open_nodes][2]=end;
	}
      }
    }
  }
  return total_var;
}

void compute_ltrees_variable_importance_multiregr(SCORE_TYPE *attribute_importance, int obj) {
  /* si obj est à -1 -> calcul classique
   * sinon, on ne parcourt que la branche de l'arbre par laquelle l'objet passe 
   * utile surtout pour un ensemble d'arbre. Permet de retrouver une sorte de branche
   * de tests (floues ou prototype) correspondant à cet objet.
   */

  int i,t,j;
  SCORE_TYPE sum_val=0.0;
  SCORE_TYPE sum_val2=0.0;
  SCORE_TYPE total_var_one_tree;
  int *ts_vector=current_learning_set;
  int length_ts_vector=global_learning_set_size;

  /*
  printf("length ls=%d\n", length_ts_vector);
  fflush(stdout);
  */

  /* remplit le vecteur */
  for (i=0; i<length_ts_vector; i++)  
    ts_vector[i]=i;

  for (i=0; i<nb_attributes; i++) {
    attribute_importance[i]=0.0;
    attribute_position[attribute_vector[i]]=i;
  }

  /* boucle sur les arbres */
  if ((current_nb_of_ensemble_terms==0) && (index_nodes>=0)) {
    /* il y a un arbre mais pas d'ensemble. On calcule l'importance de cet arbre uniquement */

    total_var_one_tree=
      compute_one_tree_variable_importance_multiregr(0, ts_vector, length_ts_vector, 1.0, attribute_importance, obj);
    sum_val2+=total_var_one_tree;
  
  } else {
    for (t=0; t<current_nb_of_ensemble_terms; t++) {

      total_var_one_tree=
	compute_one_tree_variable_importance_multiregr(ltrees[t], ts_vector, length_ts_vector, ltrees_weight[t], attribute_importance,obj);
      if (total_var_one_tree!=-1) {
	/* toutes les valeurs sont les memes */
	sum_val2=total_var_one_tree;
      }
      /* 
	 printf("%f\n", total_var_one_tree);
	 fflush(stdout);
      */
    }
  }

  /* normalizing the values */
  sum_val=0.0;
  for (i=0; i<nb_attributes; i++) {
    sum_val+=attribute_importance[i];
  }
  
  if (sum_val!=0.0) { /* arbre ne se reduit pas a la racine */
    for (i=0; i<nb_attributes; i++) {
      attribute_importance[i]/=sum_val;
    }
  }
  if (print_result) {
    printf("Total impurity(N.I(S))=%d*%f=%f\nAverage Impurity explained by trees=%f\n", 
	   length_ts_vector, 
	   sum_val2/length_ts_vector, sum_val2, 
	   sum_val/current_nb_of_ensemble_terms);
    fflush(stdout);
  }
}

/***********************************************/
/* IMPORTANCE DES VARIABLES (version multiple) */
/***********************************************/

/* code separé par rapport aux arbres classiques mais les deux pourraient être combinés sans problème */
/* pour le moment, ts_vector contient des index dans object_mapping (de 0 à N)
 */

void compute_multiregr_score_from_table_for_varimp(SCORE_TYPE *vi) {
  int i;
  SCORE_TYPE y_tot_var, n_tot_var;

  table_score_multiregr[2][0]=table_score_multiregr[0][0]-table_score_multiregr[1][0];

  for (i=0; i<nb_goal_multiregr; i++) {
    y_tot_var=0.0;
    n_tot_var=0.0;
    v_tot=0.0;

    v_tot=table_score_multiregr[0][2*i+2]-(table_score_multiregr[0][2*i+1]*table_score_multiregr[0][2*i+1])/table_score_multiregr[0][0];

    table_score_multiregr[2][2*i+1]=table_score_multiregr[0][2*i+1]-table_score_multiregr[1][2*i+1];
    table_score_multiregr[2][2*i+2]=table_score_multiregr[0][2*i+2]-table_score_multiregr[1][2*i+2];
    
    y_tot_var=fabs(table_score_multiregr[1][2*i+2]-(table_score_multiregr[1][2*i+1]*table_score_multiregr[1][2*i+1])/table_score_multiregr[1][0]);
    n_tot_var=fabs(table_score_multiregr[2][2*i+2]-(table_score_multiregr[2][2*i+1]*table_score_multiregr[2][2*i+1])/table_score_multiregr[2][0]);
    vi[i]=v_tot-(y_tot_var+n_tot_var);
  }
}

void get_vi_multiregr_separate(int *ts_vector, int start, int end, int borne, SCORE_TYPE *vi) {
  int i,j, st;
  SCORE_TYPE current_score;

  /* summarize_vector */
  for (i=0; i<2*nb_goal_multiregr+1; i++) {
    table_score_multiregr[0][i]=0.0;
    table_score_multiregr[1][i]=0.0;
  }

  for (i=start; i<=end; i++) {
    table_score_multiregr[0][0]++;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(ts_vector[i],j);
      table_score_multiregr[0][2*j+1]+=y;
      table_score_multiregr[0][2*j+2]+=y*y;
    }
  }
  
  /* calcule v_tot en fait */
  
  if ((start>=borne)||(borne>end)) {
    /* tous les objets sont a droite ou a gauche -> vi=0 
     * (ne peut arriver qu'a cause d'erreur d'arrondi ?)
     */
    for (i=0; i<nb_goal_multiregr; i++) {
      vi[i]=0;
    }
    return;
  }

  /* fill the table_score, il y a sûrement une manière plus efficace que de le faire
   * objet par objet */

  for (i=start; i<borne; i++) {
    table_score_multiregr[1][0]++;
    for (j=0; j<nb_goal_multiregr; j++) {
      float y=getobjy_multiregr_learn(ts_vector[i],j);
      table_score_multiregr[1][2*j+1]+=y;
      table_score_multiregr[1][2*j+2]+=(y*y);
    }
  }

  /* compute the score */
  compute_multiregr_score_from_table_for_varimp(vi);

}

void compute_one_tree_variable_importance_multiregr_separate(int tree, int *ts_vector, int length_ts_vector, float weight,
							     SCORE_TYPE *attribute_importance, int obj) {
  int i;
  SCORE_TYPE *vi=NULL;

  vi=(SCORE_TYPE *)malloc((size_t)nb_goal_multiregr*sizeof(SCORE_TYPE));

  if (vi==NULL) {
    printf("Impossible d'allouer de la memoire dans compute_one_tree_variable_importance_multiregr_separate\n");
    fflush(stdout);
    exit(0);
  }

  index_stack_open_nodes++;
  stack_open_nodes[index_stack_open_nodes][0]=tree;
  stack_open_nodes[index_stack_open_nodes][1]=0;
  stack_open_nodes[index_stack_open_nodes][2]=length_ts_vector-1;
  
  while (index_stack_open_nodes>=0) {
    int node=stack_open_nodes[index_stack_open_nodes][0];
    int start=stack_open_nodes[index_stack_open_nodes][1];
    int end=stack_open_nodes[index_stack_open_nodes][2];
    int node_size=end-start+1;

    if ((left_successor[node]==-1)||(node_size==1)) {
      index_stack_open_nodes--;
      /*
	printf("term start=%d, end=%d, tested attribute=%d isop=%d\n", start, end, tested_attribute[node], index_stack_open_nodes);
	fflush(stdout);
      */
    } else {
      /* separation */
      int borne=separate_ls_vector_local(tested_attribute[node], threshold[node], ts_vector, start, end);

      /*
      printf("current_node=%d, left=%d, right=%d\n",node, left_successor[node], right_successor[node]);
      printf("start=%d, borne=%d, end=%d, tested attribute=%d, threshold=%f, isop=%d\n", start, borne, end, tested_attribute[node], threshold[node].f, index_stack_open_nodes);
      fflush(stdout);
      */

      /* calcul de l'importance (seulement si borne OK) */
      get_vi_multiregr_separate(ts_vector, start, end, borne, vi);

      /*
      printf("vi[0]=%f\n",vi[0]);
      fflush(stdout);
      */

      /* mis a jour du vecteur */
      for (i=0; i<nb_goal_multiregr; i++)
	attribute_importance[i*nb_attributes+attribute_position[tested_attribute[node]]]+=(weight*vi[i]);
      
      /* left and right successors are put on the stack */
      index_stack_open_nodes--;
      if (obj<0) {
	if (start<borne) {
	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=node+left_successor[node];
	  stack_open_nodes[index_stack_open_nodes][1]=start;
	  stack_open_nodes[index_stack_open_nodes][2]=borne-1;
	}
	if (borne<=end) {
	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=node+right_successor[node];
	  stack_open_nodes[index_stack_open_nodes][1]=borne;
	  stack_open_nodes[index_stack_open_nodes][2]=end;
	}
      } else {
	if (check_test(tested_attribute[node],getattval(obj,tested_attribute[node]),threshold[node])) {
	  /* on met le gauche */
	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=node+left_successor[node];
	  stack_open_nodes[index_stack_open_nodes][1]=start;
	  stack_open_nodes[index_stack_open_nodes][2]=borne-1;
	} else {
	  /* on met le droit */
	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=node+right_successor[node];
	  stack_open_nodes[index_stack_open_nodes][1]=borne;
	  stack_open_nodes[index_stack_open_nodes][2]=end;
	}
      }
    }
  }

  free((SCORE_TYPE *)vi);

}

void compute_ltrees_variable_importance_multiregr_separate(SCORE_TYPE *attribute_importance, int obj) {
  /* si obj est à -1 -> calcul classique
   * sinon, on ne parcourt que la branche de l'arbre par laquelle l'objet passe 
   * utile surtout pour un ensemble d'arbre. Permet de retrouver une sorte de branche
   * de tests (floues ou prototype) correspondant à cet objet.
   */

  int i,t,j;
  SCORE_TYPE sum_val=0.0, sum_weight;
  SCORE_TYPE sum_val2=0.0;
  SCORE_TYPE total_var_one_tree;
  int *ts_vector=current_learning_set;
  int length_ts_vector=global_learning_set_size;

  /*
  printf("length ls=%d\n", length_ts_vector);
  fflush(stdout);
  */

  /* remplit le vecteur */
  for (i=0; i<length_ts_vector; i++)  
    ts_vector[i]=i;

  for (i=0; i<nb_attributes; i++) {
    for (j=0; j<nb_goal_multiregr; j++) {
      attribute_importance[j*nb_attributes+i]=0.0;
    }
    attribute_position[attribute_vector[i]]=i;
  }

  /* boucle sur les arbres */
  if ((current_nb_of_ensemble_terms==0) && (index_nodes>=0)) {
    /* il y a un arbre mais pas d'ensemble. On calcule l'importance de cet arbre uniquement */

    compute_one_tree_variable_importance_multiregr_separate(0, ts_vector, length_ts_vector, 1.0, attribute_importance, obj);    
  } else {
    for (t=0; t<current_nb_of_ensemble_terms; t++) {

      compute_one_tree_variable_importance_multiregr_separate(ltrees[t], ts_vector, length_ts_vector, ltrees_weight[t], attribute_importance, obj);      
      sum_weight+=ltrees_weight[t];
    }
  }

  /* normalizing the values */
  if (average_predictions_ltrees==1){
    for (i=0; i<nb_goal_multiregr*nb_attributes; i++)
      attribute_importance[i]/=sum_weight;
  }

}

/******************************/
/* TEST D'UN MODELE MULTIREGR */
/******************************/

/* définit des fonctions particulières pour l'apprentissage sur des données complexes */

/* calcul de l'erreur quadratique sur un ensemble de test */

/* retrouve l'ensemble de voisins sur un ensemble d'arbres */

/* calcul du vecteur de prédiction pour un objet (permet de calculer l'erreur quadratique, etc.
 * du côté lisp
 */

void get_multiregr_tree_pred(int tree, int obj, float *vect_pred) {
  int l=get_tree_leaf(tree, obj);
  int i, j,st,end;
  float sum_weight=0.0;

  for (j=0; j<nb_goal_multiregr; j++) {
    vect_pred[j]=0.0;
  }

  st=node_subset_desc[l][1]; end=st+node_subset_desc[l][0];

  for (i=st; i<end; i++) {
    float w=node_object_weight[i];
    sum_weight+=w;
    for (j=0; j<nb_goal_multiregr; j++) {
      vect_pred[j]+=w*getobjy_multiregr_learn(node_object_list[i],j);
    }
  }

  /* printf("pred[0]=%f, sum weight obj=%f\n",vect_pred[0],sum_weight);fflush(stdout); */
  /* on normalise */
  for (i=0; i<nb_goal_multiregr; i++) {
    vect_pred[i]/=sum_weight;
  }
}

void get_multiregr_tree_pred_savepred(int tree, int obj, float *vect_pred) {
  float *tree_vect_pred;
  int j;
  
  tree_vect_pred=get_tree_prediction_vector(tree, obj);
  
  for (j=0; j<nb_goal_multiregr; j++) {
    vect_pred[j]=tree_vect_pred[j];
  }
}

void get_multiregr_ens_pred(int obj, float *vect_pred) {
  int i, j;
  float *temp_vect_pred; float sum_weight=0.0,w;

  for (j=0; j<nb_goal_multiregr; j++) {
    vect_pred[j]=0.0;
  }
  
  temp_vect_pred=(float *)MyMalloc((size_t)nb_goal_multiregr*sizeof(float));

  for (i=0; i<current_nb_of_ensemble_terms; i++) {
    w=ltrees_weight[i];
    sum_weight+=w;
    get_multiregr_tree_pred(ltrees[i],obj,temp_vect_pred);
    for (j=0; j<nb_goal_multiregr; j++) {
      vect_pred[j]+=w*temp_vect_pred[j];
    }
  }

  /* printf("vect_pred[0]=%f, sum weight tree=%f\n",vect_pred[0],sum_weight);fflush(stdout); */
  /* on normalise */
  if (average_predictions_ltrees) {
    for (i=0; i<nb_goal_multiregr; i++) {
      vect_pred[i]/=sum_weight;
    }
  }

  MyFree((float *)temp_vect_pred);
}

void get_multiregr_ens_pred_savepred(int obj, float *vect_pred) {
  int i, j;
  float sum_weight=0.0,w, *tree_vect_pred;

  for (j=0; j<nb_goal_multiregr; j++) {
    vect_pred[j]=0.0;
  }
  
  for (i=0; i<current_nb_of_ensemble_terms; i++) {
    w=ltrees_weight[i];
    sum_weight+=w;
    tree_vect_pred=get_tree_prediction_vector(ltrees[i], obj);
    for (j=0; j<nb_goal_multiregr; j++) {
      vect_pred[j]+=w*tree_vect_pred[j];
    }
  }

  /* printf("vect_pred[0]=%f, sum weight tree=%f\n",vect_pred[0],sum_weight);fflush(stdout); */
  /* on normalise */
  if (average_predictions_ltrees) {
    for (i=0; i<nb_goal_multiregr; i++) {
      vect_pred[i]/=sum_weight;
    }
  }
}

/*********************/
/* MART en multiregr */
/*********************/

ERROR_TYPE *temp_output;
int length_temp_output=0;

int temp_best_first_max_nb_tests;

float getobjy_multiregr_learn_temp(int o, int g) {
  return (float)temp_output[g*global_learning_set_size+o];
}

float make_ls_vector_mart_multiregr(int tree) {
  int i,j, predc, truec, nb_error=0;

  if (tree<0) { /* initialisation */
    for (i=0; i<global_learning_set_size; i++)
      object_weight[i]=1.0;
    /* on initialise temp_output si ce n'est pas fait */
    if (length_temp_output<(global_learning_set_size*nb_goal_multiregr)) {
      if (length_temp_output!=0)
	MyFree((ERROR_TYPE *)temp_output);
      temp_output=(ERROR_TYPE *)MyMalloc((size_t)global_learning_set_size*nb_goal_multiregr*sizeof(ERROR_TYPE));
      if (!temp_output) {
	printf("impossible d'allouer de la memoire dans make_ls_vector_mart_multiregr\n");
	exit(0);
      }
      length_temp_output=global_learning_set_size;
    }
    
    /* on le remplit */
    for (i=0; i<global_learning_set_size; i++) {
      for (j=0; j<nb_goal_multiregr; j++) {
	temp_output[j*global_learning_set_size+i]=getobjy_multiregr_learn(i,j);
      }
    }

    /* on fixe la fonction de sortie */
    getobjy_multiregr_learn=getobjy_multiregr_learn_temp;
    /* pour que le premier modele ait une complexite de 1 (pas de test) */
    temp_best_first_max_nb_tests=best_first_max_nb_tests;
    best_first_max_nb_tests=0;

  } else {

    /* met à jour l'erreur sur chaque objet */
    if (current_nb_of_ensemble_terms==0) {
      best_first_max_nb_tests=temp_best_first_max_nb_tests;
      for (i=0; i<global_learning_set_size; i++) {
	float *tree_vect_pred=get_tree_prediction_vector(tree,object_mapping[i]);
	for (j=0; j<nb_goal_multiregr; j++) {
	  temp_output[j*global_learning_set_size+i]-=tree_vect_pred[j];
	}
      }
      return 1.0;
    } else {
      for (i=0; i<global_learning_set_size; i++) {
	float *tree_vect_pred=get_tree_prediction_vector(tree,object_mapping[i]);
	for (j=0; j<nb_goal_multiregr; j++) {
	  temp_output[j*global_learning_set_size+i]-=mart_mu*tree_vect_pred[j];
	}
      }
      return mart_mu;
    }
  }

  return 1.0;

}


int separate_ls_and_weight_vector_local(int best_attribute, union threshold_type best_threshold, int *ls_vector, float *weight_vector, int start, int end) {
  
  while (start!=end) {
    
    while (start!=end && (check_test(best_attribute,
				     getattval(object_mapping[ls_vector[start]],best_attribute), 
				     best_threshold))) {
      start++;
    }
    while (start!=end && !(check_test(best_attribute,
				      getattval(object_mapping[ls_vector[end]],best_attribute), 
				      best_threshold))) {
      end--;
    }
    if (start!=end) { /* on inverse les deux */
      int temp; float temp_f;

      temp=ls_vector[start];
      ls_vector[start]=ls_vector[end];
      ls_vector[end]=temp;

      temp_f=weight_vector[start];
      weight_vector[start]=weight_vector[end];
      weight_vector[end]=temp_f;
      
      start++;
    }
  }
  /* ici, on a start=end, on renvoie la borne */

  if (check_test(best_attribute,getattval(object_mapping[ls_vector[start]],best_attribute), best_threshold))
    return (start+1);
  else
    return start;
}


void compute_node_subset_current_ensemble_weighted(float *ow) {
  int t,i,j,s,size_ls;
  int pos_nol; /* position in node_object_list */

  object_weight2=ow;

  if (save_ensemble_ls && (save_ensemble_ls_nb_ls>0)) { /* si on a sauvegarde le LS, on le recupere */
                                  /* Notons que ça peut être mieux de reprendre le LS 
				   *  complet à chaque fois (p.e. pour le bagging)
				   */
    if (node_object_allocatedp && (node_object_weight!=NULL))
      free((float *)node_object_weight);
    if (node_object_allocatedp && (node_object_list!=NULL))
      free((int *)node_object_list);
    node_object_allocatedp=0;

    node_object_list=save_ensemble_ls_vector;
    node_object_weight=save_ensemble_ls_weight;

    /* affiche les LS */
    
    /*
    pos_nol=0;
    printf("save_ensemble_ls_pos=%d\n",save_ensemble_ls_pos);
    for (i=0;i<current_nb_of_ensemble_terms; i++) {
      for (j=0; j<save_ensemble_ls_size[i]; j++) {
	printf("%d,o=%d,w=%f ", pos_nol, node_object_list[pos_nol], save_ensemble_ls_weight[pos_nol]);
	pos_nol++;
      }
      printf("\n");
    }
    fflush(stdout);
    */
    
  } else { /* sinon, on cree tout */
    if (node_object_allocatedp && (node_object_weight!=NULL))
      free((float *)node_object_weight);
    if (node_object_allocatedp && (node_object_list!=NULL))
      free((int *)node_object_list);
    s=current_nb_of_ensemble_terms*global_learning_set_size;
    node_object_list=(int *)malloc((size_t)s*sizeof(int));
    node_object_weight=(float *)malloc((size_t)s*sizeof(float));
    node_object_allocatedp=1;

    pos_nol=0;
    for (t=0; t<current_nb_of_ensemble_terms; t++) {
      for (i=0; i<global_learning_set_size; i++) {
	node_object_list[pos_nol]=i;
	node_object_weight[pos_nol]=1.0;
	pos_nol++;
      }
    }
  }

  /*
  printf("pos_nol=%d, global_learning_set_size=%d, save_ensemble_ls=%d\n",pos_nol,global_learning_set_size,save_ensemble_ls);
  fflush(stdout);
  */

  /* alloue la mémoire pour la table node_subset_desc */
  allocate_node_subset_desc(index_nodes+1);

  /* we fill node_object_list with all objects from the ls (index to object_mapping) */
  pos_nol=0;

  for (t=0; t<current_nb_of_ensemble_terms; t++) {
        
    if (save_ensemble_ls && (save_ensemble_ls_nb_ls>0))
      size_ls=save_ensemble_ls_size[t];
    else
      size_ls=global_learning_set_size;

    index_stack_open_nodes=-1;
    index_stack_open_nodes++;
    stack_open_nodes[index_stack_open_nodes][0]=ltrees[t];
    stack_open_nodes[index_stack_open_nodes][1]=pos_nol;
    stack_open_nodes[index_stack_open_nodes][2]=pos_nol+size_ls-1;

    /* on remplit le vector */
    node_subset_desc[ltrees[t]][0]=size_ls;
    node_subset_desc[ltrees[t]][1]=pos_nol;
    
    /* 
    for (i=0; i<global_learning_set_size; i++)
      node_object_list[pos_nol++]=i;
    */
    
    pos_nol+=size_ls;
    
    while(index_stack_open_nodes>=0) {
      int node=stack_open_nodes[index_stack_open_nodes][0];
      int start=stack_open_nodes[index_stack_open_nodes][1];
      int end=stack_open_nodes[index_stack_open_nodes][2];
      int node_size=end-start+1;

      /*
	printf("node=%d,start=%d,end=%d\n",node,start,end);fflush(stdout);
	printlsvector(node_object_list, start, end);
      */
      
      index_stack_open_nodes--;

      if (left_successor[node]!=-1) { /* not a leaf */
	int borne=separate_ls_and_weight_vector_local(tested_attribute[node], threshold[node], node_object_list, node_object_weight, start, end);
	int left=node+left_successor[node], right=node+right_successor[node];

	node_subset_desc[left][0]=borne-start;
	node_subset_desc[left][1]=start;
	node_subset_desc[right][0]=end-borne+1;
	node_subset_desc[right][1]=borne;

	index_stack_open_nodes++;
	stack_open_nodes[index_stack_open_nodes][0]=node+left_successor[node];
	stack_open_nodes[index_stack_open_nodes][1]=start;
	stack_open_nodes[index_stack_open_nodes][2]=borne-1;
	index_stack_open_nodes++;
	stack_open_nodes[index_stack_open_nodes][0]=node+right_successor[node];
	stack_open_nodes[index_stack_open_nodes][1]=borne;
	stack_open_nodes[index_stack_open_nodes][2]=end;
      } /* else {
	   ici, on pourrait normaliser les poids entre start et end pour qu'il somme à 1
	   cela eviterait de devoir normaliser avant chaque prédiction
	   } */
    }
    /* printf("."); fflush(stdout); */
  }
}
