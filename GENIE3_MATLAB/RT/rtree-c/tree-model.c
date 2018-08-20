
/*
 * Author: Pierre Geurts (p.geurts@ulg.ac.be)
 * (c) 2002-2010
 *
 * (Please do not distribute without permission of the author)
 *
 */


void init_threshold_type(union threshold_type *u) {
  int i;
  for (i=0;i<MAX_NUMBER_OF_SYMBOLIC_VALUES_DIV_32; i++)
    (*u).i[i]=0;
}

void add1_threshold_type(union threshold_type *u) {
  int i=0;
  do {
    (*u).i[i]++;
  } while (((*u).i[i++]==0)&&(i<=MAX_NUMBER_OF_SYMBOLIC_VALUES_DIV_32));
}

/**********************/
/* VARIABLES GLOBALES */
/**********************/

/* DEFINITION DE LA "TABLE" POUR STOCKER LES ARBRES */

int size_current_tree_table=0;
int size_current_tree_table_pred=0;

union threshold_type *threshold=NULL; /* le seuil associe a un noeud (sur 64 BITS pour représenter plus de valeurs
				                      * symboliques */
int *left_successor=NULL; /* la position RELATIVE du successeur par rapport */
int *right_successor=NULL; /* a ce noeud */
int *tested_attribute=NULL; /* l'attribut teste en ce noeud */
int *prediction=NULL; /* l'index de la prediction associe au noeud */
float *node_weight=NULL; /* un poids pour le test en dual perturb and combine */
/* Pour missing values: la taille de la feuille en nb d'elements */
float *node_size=NULL;

float **prediction_values=NULL;
int index_prediction_values=-1;

int index_nodes=-1;

int stack_open_nodes[MAX_OPEN_NODES][3];

/* pour le best first */
SCORE_TYPE priority_open_nodes[MAX_OPEN_NODES];
/* pour le true best first */
int *tbf_best_attribute=NULL;
int *tbf_stop_splitting=NULL;
int *tbf_not_significant=NULL;
SCORE_TYPE *tbf_best_threshold_score=NULL;
SCORE_TYPE *tbf_best_threshold_info=NULL;
union threshold_type *tbf_best_threshold=NULL;
int size_tbf=0;

int index_stack_open_nodes=-1;

/* table des valeurs */

CORETABLE_TYPE *core_table=NULL;

int nb_obj_in_core_table;
/* vector with attribute descriptor (0 if numerical, the number of values
 * if symbolic)
 */
int *attribute_descriptors=NULL;
int length_attribute_descriptors;

/* parametres generaux pour la construction */

int goal; /* goal_classification ou goal_regression */

int *attribute_vector=NULL; /* indice des attributs */
int nb_attributes;

int *current_learning_set=NULL; /* l'ensemble d'apprentissage */
int current_learning_set_size;
int *object_mapping=NULL;
SCORE_TYPE *object_weight=NULL;
int global_learning_set_size=0; /* taille des table object_mapping et object_weight: peut être plus grand
				   que current_learning_set_size si certains objets sont de poids nuls*/

int min_node_size; /* nombre min d'objets dans un noeud */

int best_first=0; /* si 1 arbre construit en best_first (sur base de la variance dans le noeud) */
int true_best_first=0; /* si arbre construit en vrai best first (sur base du score du meilleur test) */
                       /* on suppose que true_best_first a priorité sur best_first si les deux sont vrais */

int best_first_max_nb_tests=1; /* nb maximum de tests si en (true) best first */

/* parametres du tri: *
 * par defaut, on trie chaque fois localement */
void (*sort_ls_vector)(int *ls_vector, int start, int end, int att)=quicksort_ls_vector;
int (*separate_ls_vector)(int best_attribute, union threshold_type best_threshold, int *ls_vector, int start, int end)=separate_ls_vector_local; 

/* methodes d'ensemble */

float (*make_ls_vector)(int tree);
int number_of_ensemble_terms=50;
int save_ensemble_while_growing=0;
int store_ensemble=1;

/* Niveau d'affichage de messages a l'ecran (0 -> rien, 1 -> tout) */
int print_result=1;

/* parametres generaux pour la classification */

int nb_classes; /* nombre de classes pour la classification */

/* parametres generaux pour la regression */ 

/* ... */

/* pour le calcul du score */

SCORE_TYPE **table_score=NULL; /* [3][MAX_NUMBER_OF_PREDICTION_VALUES+1]; */
SCORE_TYPE **table_score_symb=NULL; /* [MAX_NUMBER_OF_SYMBOLIC_VALUES][MAX_NUMBER_OF_PREDICTION_VALUES+1]; */
int nb_of_predictions;

/* on pourrait utiliser table_score et table_score_symb mais je préfère utiliser
 * des tables dédiées pour que MAX_GOAL_MULTIREGR soit indépendant de MAX_NUMBER_OF_NODES
 */
SCORE_TYPE **table_score_multiregr=NULL; /* [3][2*MAX_GOAL_MULTIREGR+1]; */
SCORE_TYPE **table_score_symb_multiregr=NULL; /* [MAX_NUMBER_OF_SYMBOLIC_VALUES][2*MAX_GOAL_MULTIREGR+1]; */

/* parametre de stop splitting pour les arbres */

float h_min;

/* variables globales pour le calcul des meilleures score */

int best_attribute;
union threshold_type best_threshold;
SCORE_TYPE best_threshold_score;

union threshold_type current_threshold;
SCORE_TYPE current_threshold_score;

/* variables globales pour le stop splitting criterio a posteriori (not-significant) */
SCORE_TYPE current_threshold_info;
SCORE_TYPE best_threshold_info;

/* pour le calcul du score */
SCORE_TYPE info;
SCORE_TYPE v_tot;
SCORE_TYPE v_min=0.0;
float alfa_rt=1.0;


/***************/
/* CODE ARBRES */
/***************/

/* BUILD_ONE_TREE */
/*----------------*/
/* fonction generale de construction d'un arbre */

/* print_stack */
void print_stack(int size) {
  int i;
  printf("Begin stack:\n");
  for (i=0; i<=size; i++) {
    printf("Pos %d - node=%d (%d-%d) - score=%e - att=%d\n",i,stack_open_nodes[i][0], stack_open_nodes[i][1], stack_open_nodes[i][2], tbf_best_threshold_info[i], tbf_best_attribute[i]);
  }
  printf("End stack\n");
  fflush(stdout);
}

int build_one_tree() {
  int tree, i;
  int nb_tests=0;

  /* on construit un noeud et on le place sur le stack */
  index_nodes++;
  prediction[index_nodes]=-1; /* par defaut, c'est une feuille */
  tested_attribute[index_nodes]=-1; /* pas de test */
  left_successor[index_nodes]=-1;
  right_successor[index_nodes]=-1;
  tree=index_nodes;
  
  index_stack_open_nodes++;
  stack_open_nodes[index_stack_open_nodes][0]=tree;
  stack_open_nodes[index_stack_open_nodes][1]=0;
  stack_open_nodes[index_stack_open_nodes][2]=current_learning_set_size-1;
  
  if (true_best_first) {
    /* on recherche le meilleur split et on stocke l'info */
    summarize_vector(current_learning_set, 0, current_learning_set_size-1);

    node_size[tree]=table_score[0][0]; /* valeurs manquantes ? */
    
    tbf_stop_splitting[index_stack_open_nodes]=((current_learning_set_size==1) || (current_learning_set_size<min_node_size) || stop_splitting_criterio());
    if (!tbf_stop_splitting[index_stack_open_nodes]) {
      find_a_split(current_learning_set, 0, current_learning_set_size-1);
      tbf_best_attribute[index_stack_open_nodes]=best_attribute;
      /* tbf_best_threshold_score[index_stack_open_nodes]=best_threshold_score; */
      tbf_best_threshold_info[index_stack_open_nodes]=best_threshold_info;
      tbf_best_threshold[index_stack_open_nodes]=best_threshold;
      tbf_not_significant[index_stack_open_nodes]=not_significant_test();
    } else {
      tbf_best_threshold_info[index_stack_open_nodes]=-1.0;
    }
  } else if (best_first) {
    /* on calcule la variance totale */
    summarize_vector(current_learning_set, 0, current_learning_set_size-1);
    priority_open_nodes[index_stack_open_nodes]=best_first_node_priority_function();
  }

  /* on lance la boucle de developpement, tant qu'il y a des noeuds ouverts */
  while (index_stack_open_nodes>=0) {
    int node=stack_open_nodes[index_stack_open_nodes][0];
    int start=stack_open_nodes[index_stack_open_nodes][1];
    int end=stack_open_nodes[index_stack_open_nodes][2];
    int nodesize=end-start+1;
    
    if (!true_best_first) {
      /* resume du vecteur */
      summarize_vector(current_learning_set, start, end);
      /* pour missing_values */
      node_size[node]=table_score_multiregr[0][0];
    }

    /* condition d'arret: */
    if ((true_best_first && (tbf_stop_splitting[index_stack_open_nodes] || (nb_tests>=best_first_max_nb_tests)))||
	(!true_best_first && ((best_first && (nb_tests>=best_first_max_nb_tests)) ||
			      (nodesize==1) ||
			      (nodesize<min_node_size) || 
			      stop_splitting_criterio()))) {
      /* c'est une feuille */
      if (true_best_first) {
	summarize_vector(current_learning_set, start, end);
      }
      prediction[node]=make_leaf_prediction();
      index_stack_open_nodes--;

    } else {

      if (!true_best_first) {
	/* on recherche le meilleur split */
	find_a_split(current_learning_set, start, end);
      } else {
	/* on reprend l'info qu'on avait stockée */
	best_attribute=tbf_best_attribute[index_stack_open_nodes];
	best_threshold=tbf_best_threshold[index_stack_open_nodes];
	best_threshold_score=tbf_best_threshold_score[index_stack_open_nodes];
	best_threshold_info=tbf_best_threshold_info[index_stack_open_nodes];
      }
	
      /* le resultat est dans best_threshold_score), best_attribute */
      if ((true_best_first && tbf_not_significant[index_stack_open_nodes])
	  ||(!true_best_first && not_significant_test())) {
	/* c'est une feuille */
	/* printf("feuille");fflush(stdout); */
	prediction[node]=make_leaf_prediction();
	index_stack_open_nodes--;
      } else {
	/* on separe les objets sur base */
	int left,right,borne;

	borne=separate_ls_vector(best_attribute, best_threshold, current_learning_set,start,end);

	if ((borne<=start) || (borne>end)) {
 	
	  printf("ZERO\nstart=%d, borne=%d, end=%d\n", start, borne, end);fflush(stdout);
	  printf("best_attribute=%d, best_threshold=%f, best_threshold_score=%e best_info=%e\n",best_attribute, best_threshold.f,
		 best_threshold_score, best_threshold_info);
	  fflush(stdout);	
	}

	nb_tests++;

	/* on cree deux nouveaux noeuds */
	index_nodes++; left=index_nodes;
	index_nodes++; right=index_nodes;
	prediction[left]=-1; prediction[right]=-1;
	tested_attribute[left]=-1; tested_attribute[right]=-1;
	left_successor[left]=-1; left_successor[right]=-1;
	right_successor[left]=-1; right_successor[right]=-1;

	/* on met a jour le noeud courant */
	threshold[node]=best_threshold;
	tested_attribute[node]=best_attribute;
	left_successor[node]=left-node;
	right_successor[node]=right-node;
	
	/* on place les nouveaux noeuds sur la pile */
	/* La manière de les placer depend de true_best_first et best_first */

	if (true_best_first) { /* on trie les noeuds par ordre décroissant du réduction de variance (ou info) */
	  int i; int ns; int stsp;
	  float pr_val;

	  /* left successor */
	  ns=borne-start;
	  summarize_vector(current_learning_set, start, borne-1);

	  node_size[node]=table_score_multiregr[0][0];

	  stsp=((ns==1) || (ns<min_node_size) || stop_splitting_criterio());
	  if (!stsp) {
	    find_a_split(current_learning_set, start, borne-1);	  
	    pr_val=best_threshold_info;
	  } else {
	    pr_val=-1.0;
	  }
	  
	  i=index_stack_open_nodes-1;
	  while ((i>=0) && (tbf_best_threshold_info[i]>pr_val)) {
	      /* on pousse le noeud jusqu'à sa place */
	      stack_open_nodes[i+1][0]=stack_open_nodes[i][0];
	      stack_open_nodes[i+1][1]=stack_open_nodes[i][1];
	      stack_open_nodes[i+1][2]=stack_open_nodes[i][2];

	      tbf_best_attribute[i+1]=tbf_best_attribute[i];
	      tbf_best_threshold_score[i+1]=tbf_best_threshold_score[i];
	      tbf_best_threshold_info[i+1]=tbf_best_threshold_info[i];
	      tbf_best_threshold[i+1]=tbf_best_threshold[i];
	      tbf_not_significant[i+1]=tbf_not_significant[i];
	      tbf_stop_splitting[i+1]=tbf_stop_splitting[i];

	      i--;
	  } 
	  stack_open_nodes[i+1][0]=left;
	  stack_open_nodes[i+1][1]=start;
	  stack_open_nodes[i+1][2]=borne-1;

	  tbf_stop_splitting[i+1]=stsp;
	  tbf_best_threshold_info[i+1]=pr_val;
	  if (!stsp) {
	    tbf_best_attribute[i+1]=best_attribute;
	    tbf_best_threshold_score[i+1]=best_threshold_score;
	    tbf_best_threshold[i+1]=best_threshold;
	    tbf_not_significant[i+1]=not_significant_test();
	  }

	  /* right successor */
	  index_stack_open_nodes++;
	  ns=end-borne+1;
	  summarize_vector(current_learning_set, borne, end);

	  node_size[node]=table_score_multiregr[0][0];

	  stsp=((ns==1) || (ns<min_node_size) || stop_splitting_criterio());
	  if (!stsp) {
	    find_a_split(current_learning_set, borne, end);	  
	    pr_val=best_threshold_info;
	  } else {
	    pr_val=-1.0;
	  }
	  
	  i=index_stack_open_nodes-1;
	  while ((i>=0) && (tbf_best_threshold_info[i]>pr_val)) {
	      /* on pousse le noeud jusqu'à sa place */
	      stack_open_nodes[i+1][0]=stack_open_nodes[i][0];
	      stack_open_nodes[i+1][1]=stack_open_nodes[i][1];
	      stack_open_nodes[i+1][2]=stack_open_nodes[i][2];

	      tbf_best_attribute[i+1]=tbf_best_attribute[i];
	      tbf_best_threshold_score[i+1]=tbf_best_threshold_score[i];
	      tbf_best_threshold_info[i+1]=tbf_best_threshold_info[i];
	      tbf_best_threshold[i+1]=tbf_best_threshold[i];
	      tbf_not_significant[i+1]=tbf_not_significant[i];
	      tbf_stop_splitting[i+1]=tbf_stop_splitting[i];
	      
	      i--;
	  } 
	  stack_open_nodes[i+1][0]=right;
	  stack_open_nodes[i+1][1]=borne;
	  stack_open_nodes[i+1][2]=end;
	  
	  tbf_stop_splitting[i+1]=stsp;
	  tbf_best_threshold_info[i+1]=pr_val;
	  if (!stsp) {
	    tbf_best_attribute[i+1]=best_attribute;
	    tbf_best_threshold_score[i+1]=best_threshold_score;
	    tbf_best_threshold[i+1]=best_threshold;
	    tbf_not_significant[i+1]=not_significant_test();
	  }

	} else if (best_first) { /* on trie les noeuds par ordre de priorite. */
	  int i;
	  float pr_val;

	  /* left successor */
	  summarize_vector(current_learning_set, start, borne-1);
	  pr_val=best_first_node_priority_function();
	  i=index_stack_open_nodes-1;
	  while ((i>=0) && (priority_open_nodes[i]>pr_val)) {
	      /* on pousse le noeud jusqu'à sa place */
	      stack_open_nodes[i+1][0]=stack_open_nodes[i][0];
	      stack_open_nodes[i+1][1]=stack_open_nodes[i][1];
	      stack_open_nodes[i+1][2]=stack_open_nodes[i][2];
	      priority_open_nodes[i+1]=priority_open_nodes[i];
	      i--;
	  } 
	  stack_open_nodes[i+1][0]=left;
	  stack_open_nodes[i+1][1]=start;
	  stack_open_nodes[i+1][2]=borne-1;
	  priority_open_nodes[i+1]=pr_val;

	  /* right successor */
	  index_stack_open_nodes++;
	  summarize_vector(current_learning_set, borne, end);
	  pr_val=best_first_node_priority_function();
	  i=index_stack_open_nodes-1;
	  while ((i>=0) && (priority_open_nodes[i]>pr_val)) {
	      /* on pousse le noeud jusqu'à sa place */
	      stack_open_nodes[i+1][0]=stack_open_nodes[i][0];
	      stack_open_nodes[i+1][1]=stack_open_nodes[i][1];
	      stack_open_nodes[i+1][2]=stack_open_nodes[i][2];
	      priority_open_nodes[i+1]=priority_open_nodes[i];
	      i--;
	  } 
	  stack_open_nodes[i+1][0]=right;
	  stack_open_nodes[i+1][1]=borne;
	  stack_open_nodes[i+1][2]=end;
	  priority_open_nodes[i+1]=pr_val;


	} else { /* pas de best_first, on les met betement sur la pile */

	  stack_open_nodes[index_stack_open_nodes][0]=left;
	  stack_open_nodes[index_stack_open_nodes][1]=start;
	  stack_open_nodes[index_stack_open_nodes][2]=borne-1;
	  index_stack_open_nodes++;
	  stack_open_nodes[index_stack_open_nodes][0]=right;
	  stack_open_nodes[index_stack_open_nodes][1]=borne;
	  stack_open_nodes[index_stack_open_nodes][2]=end;

	}	  
      }     
    }
  }
  return tree;
}

int separate_ls_vector_local(int best_attribute, union threshold_type best_threshold, int *ls_vector, int start, int end) {
  
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
      int temp;
      temp=ls_vector[start];
      ls_vector[start]=ls_vector[end];
      ls_vector[end]=temp;
      start++;
    }
  }
  /* ici, on a start=end, on renvoie la borne */

  if (check_test(best_attribute,getattval(object_mapping[ls_vector[start]],best_attribute), best_threshold))
    return (start+1);
  else
    return start;
}

/* fonction de tri (from the numerical recipes in C) */

#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp
#define M_QS 7
#define QUICK_SORT_STACK_SIZE 50
#define VAL(a) getattval(object_mapping[a],att)

/* Here M is the size of subarrays sorted by straight insertion and NSTACK is the required auxiliary storage. */

unsigned long *istack[QUICK_SORT_STACK_SIZE];

void quicksort_ls_vector(int *ls_vector, int start, int end, int att) {

  /* Sorts an array arr[1..n] into ascending numerical order using the Quicksort algorithm. n is input; 
     arr is replaced on output by its sorted rearrangement. (extrait de numerical recipes in c) */

  int i,ir=end,j,k,l=start;
  int istack[QUICK_SORT_STACK_SIZE];
  int jstack=-1,o;
  int temp;
  float a;
  
  for (;;) {
    /* Insertion sort when subarray small enough.*/ 
    if (ir-l < M_QS) {
      for (j=l+1;j<=ir;j++) {
	o=ls_vector[j];
	a=VAL(o);
	for (i=j-1;i>=l;i--) {
	  if (VAL(ls_vector[i]) <= a)
	    break;
	  ls_vector[i+1]=ls_vector[i];
	}
	ls_vector[i+1]=o; 
      } 
      if (jstack == -1) 
	break; 
      ir=istack[jstack--];  /* Pop stack and begin a new round of partitioning. */ 
      l=istack[jstack--]; 
    } else { 
      k=(l+ir) >> 1; /* Choose median of left, center, and right elements as partitioning element a. */
                     /* Also rearrange so that a[l]<=a[l+1]<=a[ir]. */

      SWAP(ls_vector[k],ls_vector[l+1]);
      if (VAL(ls_vector[l]) > VAL(ls_vector[ir])) {
	SWAP(ls_vector[l],ls_vector[ir]);
      }
      if (VAL(ls_vector[l+1]) > VAL(ls_vector[ir])) {
	SWAP(ls_vector[l+1],ls_vector[ir]);
      }
      if (VAL(ls_vector[l]) > VAL(ls_vector[l+1])) {
	SWAP(ls_vector[l],ls_vector[l+1]);
      }
      i=l+1; /* Initialize pointers for partitioning. */
      j=ir;
      o=ls_vector[l+1];
      a=VAL(o); /* Partitioning element. */ 
      for (;;) { /* Beginning of innermost loop.*/ 
	do i++; while (VAL(ls_vector[i]) < a); /* Scan up to  nd element > a. */ 
	do j--; while (VAL(ls_vector[j]) > a); /* Scan down to  nd element < a. */ 
	if (j < i) 
	  break; /* Pointers crossed. Partitioning complete. */
	SWAP(ls_vector[i],ls_vector[j]); /* Exchange elements. */ 
      } /* End of innermost loop. */
      ls_vector[l+1]=ls_vector[j]; /* Insert partitioning element. */
      ls_vector[j]=o;
      jstack += 2; /* Push pointers to larger subarray on stack, process smaller subarray immediately. */
      if (jstack > QUICK_SORT_STACK_SIZE) {
	printf("Stack too small in quicksort.\n");
	return;
      }
      if (ir-i+1 >= j-l) { 
	istack[jstack]=ir; istack[jstack-1]=i; ir=j-1; 
      } else { 
	istack[jstack]=j-1; istack[jstack-1]=l; l=i;
      }
    }
  }
}

/* DISPATCH ON ATTRIBUTE TYPE */
/*----------------------------*/

void find_a_threshold(int att, int *ls_vector, int start, int end) {
  /* simple dispatch sur le type de l'attribut */
  if (NUMERICALp(att)) {
    find_a_threshold_num(att, ls_vector, start, end);
  } else if (SYMBOLICp(att)) {
    find_a_threshold_symb(att, ls_vector, start, end);
  }
}

/* CHECK A TEST */
/* ------------ */

int check_test(int att, float val, union threshold_type threshold) {
  if (NUMERICALp(att)) {
    return (val<threshold.f);
  } else {
    return (BITN(threshold,(int)val));
  }
}

/* ARBRES CLASSIQUES */
/* ----------------- */

/* GENERAL */

/* find_the_best_split */
/* recherche le meilleur split sur tous les attributs */

int print_detailled_result=0;

void find_the_best_split(int *ls_vector, int start, int end) {
  int i;

  best_attribute=-1;
  best_threshold_score=-1.0;
  best_threshold_info=-1.0;

  if (print_detailled_result) {
    printf("Node %d (s=%d -> e=%d)\n", stack_open_nodes[index_stack_open_nodes][0], start, end);
    fflush(stdout);
  }

  for(i=0; i<nb_attributes; i++) { /* on considere tous les attributs */

    if (print_detailled_result) {
      if ((i%100)==0) {
	printf("%d...",i);
	fflush(stdout);
      }
    }
    find_a_threshold(attribute_vector[i], ls_vector, start, end);

    if (print_detailled_result) {
      printf("att %d, threshold %f, score %e\n", attribute_vector[i], current_threshold.f, current_threshold_score);
      fflush(stdout);
    }
    
    if ((current_threshold_score>=0.0) && (current_threshold_score>best_threshold_score)) {
      best_threshold_score=current_threshold_score;
      best_threshold_info=current_threshold_info;
      best_threshold=current_threshold;
      best_attribute=attribute_vector[i];
    }
  }
  
  if (print_detailled_result) {
    printf("\n");
    printf("best split: Att %d - Score %e \n", best_attribute, best_threshold_score); 
    fflush(stdout);
    fflush(stdout);
  }

}

/***************
 * EXTRA-TREES *
 ***************/

float random_split_score_threshold=0.1;
int nb_of_random_tests=1;

void find_a_split_at_random_et(int *ls_vector, int start, int end) {
  int random_att_pos, temp, nb_try=0;
  int remaining_att=nb_attributes;

  best_attribute=-1;
  best_threshold_score=-1.0;
  best_threshold_info=-1.0;

  if (print_detailled_result) {
    printf("Node %d (s=%d -> e=%d)\n", stack_open_nodes[index_stack_open_nodes][0], start, end);
    fflush(stdout);
  }
  
  nb_try=0;
  do {
    nb_try++;
    random_att_pos=get_random_integer(remaining_att);
    if (random_att_pos>=remaining_att) {
      printf("erreur dans le generateur de nombres aleatoires\n");
    }

    find_a_threshold(attribute_vector[random_att_pos], ls_vector, start, end);

    if (print_detailled_result) {
      printf("att %d, threshold %f, score %e\n", attribute_vector[random_att_pos], current_threshold.f, current_threshold_score);
      fflush(stdout);
    }

    if (current_threshold_score>best_threshold_score) {
      best_threshold_score=current_threshold_score;
      best_threshold_info=current_threshold_info;
      best_threshold=current_threshold;
      best_attribute=attribute_vector[random_att_pos];
    }

    remaining_att--;
    if (remaining_att!=0) {
      temp=attribute_vector[random_att_pos];
      attribute_vector[random_att_pos]=attribute_vector[remaining_att];
      attribute_vector[remaining_att]=temp;
    }
    
    if (current_threshold_score<0.0)
      /* l'attribut etait constant, ca ne compte pas */
      nb_try--; 
    
  } while ((remaining_att!=0) && (best_threshold_score<random_split_score_threshold) && (nb_try<nb_of_random_tests));

  if (print_detailled_result) {
    printf("\n");
    printf("best split: Att %d - Score %e \n", best_attribute, best_threshold_score); 
    fflush(stdout);
    fflush(stdout);
  }
}

/********************
 * RANDOM FORESTS   *
 ********************/

/* implementation de la variante de Breiman pour les arbres aleatoires
 * on selectionne un certain nombre d'attributs aleatoirement dans l'ensemble
 * des attributs candidats pour lesquels on recherche les meilleurs splits
 */

int rf_k=1;

void find_the_best_split_among_k(int *ls_vector, int start, int end) {
  int i;
  int remaining_att=nb_attributes;
  int random_att_pos, temp;

  best_attribute=-1;
  best_threshold_score=-1.0;
  best_threshold_info=-1.0;

  for (i=0; (i<rf_k)&&(remaining_att!=0) ; i++) {
    random_att_pos=get_random_integer(remaining_att);
    if (random_att_pos>=remaining_att) {
      printf("erreur dans le generateur de nombres aleatoires\n");
    }

    find_a_threshold(attribute_vector[random_att_pos], ls_vector, start, end);

    if (current_threshold_score>best_threshold_score) {
      best_threshold_score=current_threshold_score;
      best_threshold_info=current_threshold_info;
      best_threshold=current_threshold;
      best_attribute=attribute_vector[random_att_pos];
    }
    
    /* swap les attributs dans la liste */
    remaining_att--;
    if (remaining_att!=0) {
      temp=attribute_vector[random_att_pos];
      attribute_vector[random_att_pos]=attribute_vector[remaining_att];
      attribute_vector[remaining_att]=temp;
    }
  }
}


/******************************************************************************/
/* FONCTIONS UTILITAIRES DIVERSES */

/* valable uniquement si max_val est faible par rapport a rand_max */
/* sinon, il vaut mieux utiliser get_random_long_long */

int get_random_integer(int max_val) {
  return (int)floor((double)rand()*max_val*1.0/(RAND_MAX+1.0));
}

float get_random_float() {
  return (float)((double)rand()*1.0/(RAND_MAX+1.0));
}

void get_random_subset_non_empty_full(union threshold_type *rs, int nb_val) {
  int nb_try=0,i,  nb_1;

  /* printf("%Lu\n",subsetmax); */

  if (nb_val<5) {
    init_threshold_type(rs);
    (*rs).i[0]=1+get_random_integer((1<<nb_val)-2);
  } else {
    do {
      init_threshold_type(rs);
      nb_1=0;
      for (i=0;i<nb_val; i++) {
	if (get_random_float()<0.5) {
	  SET_BITN((*rs),i);
	  nb_1++;
	}
      }
      nb_try++;
    } while (((nb_1==0)||(nb_1==nb_val))&&(nb_try<10));
  }
}


float get_random_normal() {
  float x1=get_random_float();
  float x2=get_random_float();
  return (sqrt(-2.0*log(x1))*cos(2*VALUE_PI*x2));
}

const float chi_table[12][6]={{0.99995, 16.5, 20.0, 22.7, 25.0, 27.5},
			   {0.9999, 15.2, 18.5, 21.2, 23.6, 25.8},
			   {0.99975, 13.4, 16.6, 19.2, 21.6, 23.8},
			   {0.9995, 12.1, 15.2, 17.7, 20.0, 22.1},
			   {0.999, 10.8, 13.8, 16.3, 18.5, 20.5}, 
			   {0.995, 7.88, 10.6, 12.8, 14.9, 16.7},
			   {0.99, 6.63, 9.21, 11.3, 13.3, 15.1},
			   {0.975, 5.02, 7.38, 9.35, 11.1, 12.8},
			   {0.95, 3.84, 5.99, 7.81, 9.49, 11.1},
			   {0.9, 2.71, 4.61, 6.25, 7.78, 9.24},
			   {0.75, 1.32, 2.77, 4.11, 5.39, 6.63},
			   {0.00, .000000001, .000000001, .000000001, .000000001, .000000001}};

float chi_square(int ddl, float conf) {
  int i=0;
  /* recherche de la ligne */
  while ((i<12) && (fabs(chi_table[i][0]-conf)>0.00000000001))
    i++;
  if (i==12)/* si la valeur n'est pas bonne, tout est significatif */
    i=11;
  if (ddl<=5)
    return chi_table[i][ddl];
  else
    return chi_table[i][5];
}

/* accede a la table de valeurs */

float getattval(int obj, int att) {
  return (float)core_table[att*nb_obj_in_core_table+obj];
}

void copy_matrix_a_to_b(int *a, int *b, int size) {
  int i;
  for (i=0; i<size; i++)
    b[i]=a[i];
}

void print_bit_repr(int a) {
  int i;
  int mask=1;
  for (i=0; i<32; i++) {
    printf("%d",a&mask? 1 : 0);
    fflush(stdout);
    mask=mask<<1;
  }
}

void print_n_spaces(int n) {
  int i;
  for (i=0; i<n; i++) {
    printf(" ");
  }
  fflush(stdout);
}

/* fonction de randomization d'un vecteur */
/* melange les elements d'un vecteur (pour la validation croisee) */

void randomize_vector(int *vector, int length_vector) {
  int i,temp;
  for (i=0; i<length_vector; i++) {
    int rn=get_random_integer(length_vector-i);
    temp=vector[i];
    vector[i]=vector[i+rn];
    vector[i+rn]=temp;
  }
}

/* UTILISATION DES ARBRES */

/*********************
 * PRETTY-PRINT TREE *
 *********************/

void prettyprinttree(int tree, int max_depth) {
  int current_depth, current_node, i;
  
  index_stack_open_nodes=-1;
  index_stack_open_nodes++;
  stack_open_nodes[index_stack_open_nodes][0]=tree;
  stack_open_nodes[index_stack_open_nodes][1]=0;
  
  while(index_stack_open_nodes>=0) {
    current_node=stack_open_nodes[index_stack_open_nodes][0];
    current_depth=stack_open_nodes[index_stack_open_nodes][1];
    index_stack_open_nodes--;
    if (current_node==-1) { /* else */
      print_n_spaces(current_depth);
      printf("} else {\n");
      fflush(stdout);
    } else if (current_node==-2) { /* close */
      print_n_spaces(current_depth);
      printf("}\n");
      fflush(stdout);
    } else { /* noeud normal */
      
      if (left_successor[current_node]==-1) { /* noeud terminal */
	print_n_spaces(current_depth);
	printf("return (");
	for (i=0; i<nb_classes; i++) {
	  printf("%f", prediction_values[prediction[current_node]][i]);
	  if (i==nb_classes-1)
	    printf(");\n");
	  else
	    printf(",");
	}
	fflush(stdout);
      } else if (current_depth==max_depth) { /* trop bas */
	print_n_spaces(current_depth);
	printf("...\n");
	fflush(stdout);
      } else { /* noeud test */
	print_n_spaces(current_depth);
	if (NUMERICALp(tested_attribute[current_node])) {
	  printf("if (ATT%d < %f) {\n", tested_attribute[current_node], threshold[current_node].f);
	} else {
	  printf("if (ATT%d in (", tested_attribute[current_node]);
	  for (i=0; i<attribute_descriptors[tested_attribute[current_node]]; i++) {
	    if (BITN(threshold[current_node],i)) {
	      printf("%d,",i);
	    }
	  }
	  printf(")) {\n");
	  fflush(stdout);
	}
	index_stack_open_nodes++;
	stack_open_nodes[index_stack_open_nodes][0]=-2;
	stack_open_nodes[index_stack_open_nodes][1]=current_depth;
	index_stack_open_nodes++;
	stack_open_nodes[index_stack_open_nodes][0]=current_node+right_successor[current_node];
	stack_open_nodes[index_stack_open_nodes][1]=current_depth+1;
	index_stack_open_nodes++;
	stack_open_nodes[index_stack_open_nodes][0]=-1;
	stack_open_nodes[index_stack_open_nodes][1]=current_depth;
	index_stack_open_nodes++;
	stack_open_nodes[index_stack_open_nodes][0]=current_node+left_successor[current_node];	
	stack_open_nodes[index_stack_open_nodes][1]=current_depth+1;
      }
    }
  }
}

/****************************************
 * propagation d'un objet dans un arbre *
 ****************************************/

/* application d'un arbre a un objet */

float *get_tree_prediction_vector_classical(int tree, int obj) {
  int current_node=tree;
  while (left_successor[current_node]!=-1) {
    if (check_test(tested_attribute[current_node],
		   getattval(obj,tested_attribute[current_node]), 
		   threshold[current_node]))
      current_node+=left_successor[current_node];
    else
      current_node+=right_successor[current_node];
  }
  return prediction_values[prediction[current_node]];
}

/* idem mais renvoit la feuille plutot que la prediction */

int get_tree_leaf(int tree, int obj) {
  int current_node=tree;
  while (left_successor[current_node]!=-1) {
    if (check_test(tested_attribute[current_node],
		   getattval(obj,tested_attribute[current_node]), 
		   threshold[current_node]))
      current_node+=left_successor[current_node];
    else
      current_node+=right_successor[current_node];
  }
  return current_node;
}

/************************************************************************
 * INTERFACE AVEC LE LISP *
 **************************/

/* efface tous les arbres */

void clean_all_trees() {
  /* efface toutes les tables */
  index_nodes=-1;
  index_prediction_values=-1;
  index_stack_open_nodes=-1;
  clean_all_ensemble();
}

/* alloue toutes les tables de données de manière dynamique */

SCORE_TYPE **allocate_table_score_type(int nbl, int nbc) {
  SCORE_TYPE **tableau;
  int i,j;

  tableau=(SCORE_TYPE **)MyMalloc(nbl*sizeof(SCORE_TYPE *));
  if (tableau==NULL)
    return NULL;
  for (i=0;i<nbl;i++) {
    tableau[i]=(SCORE_TYPE *)MyMalloc(nbc*sizeof(SCORE_TYPE));
    if (tableau[i]==NULL) {
      for (j=0; j<i; j++) {
	MyFree((SCORE_TYPE *)tableau[j]);
      }
      return NULL;
    }
  }
  return tableau;
}

void free_table_score_type(SCORE_TYPE **tableau, int nbl) {
  int i;

  if (tableau!=NULL) {
    for (i=0;i<nbl;i++) {
      if (tableau[i]!=NULL)
	MyFree((SCORE_TYPE *)tableau[i]);
    }
    MyFree((SCORE_TYPE **)tableau);
  }
}

float **allocate_table_float(int nbl, int nbc) {
  float **tableau;
  int i,j;

  tableau=(float **)MyMalloc(nbl*sizeof(float *));
  if (tableau==NULL)
    return NULL;
  for (i=0;i<nbl;i++) {
    tableau[i]=(float *)MyMalloc(nbc*sizeof(float));
    if (tableau[i]==NULL) {
      for (j=0; j<i; j++) {
	MyFree((float *)tableau[j]);
      }
      return NULL;
    }
  }
  return tableau;
}

void free_table_float(float **tableau, int nbl) {
  int i;

  if (tableau!=NULL) {
    for (i=0;i<nbl;i++) {
      if (tableau[i]!=NULL)
	MyFree((float *)tableau[i]);
    }
    MyFree((float **)tableau);
  }
}


void free_tree_tables() {
  int i;

  clean_all_trees();

  if (left_successor!=NULL) {
    MyFree((int *)left_successor);
    left_successor=NULL;
  }
  if (right_successor!=NULL) {
    MyFree((int *)right_successor);
    right_successor=NULL;
  }
  if (tested_attribute!=NULL) {
    MyFree((int *)tested_attribute);
    tested_attribute=NULL;
  }
  if (prediction!=NULL) {
    MyFree((int *)prediction);
    prediction=NULL;
  }
  if (node_weight!=NULL) {
    MyFree((float *)node_weight);
    node_weight=NULL;
  }
  if (node_size!=NULL) {
    MyFree((float *)node_size);
    node_size=NULL;
  }
  if (threshold!=NULL) {
    MyFree((union threshold_type *)threshold);
    threshold=NULL;
  }

  free_table_float(prediction_values,size_current_tree_table_pred);
  prediction_values=NULL;

  free_table_score_type(table_score,3);
  table_score=NULL;

  free_table_score_type(table_score_symb,MAX_NUMBER_OF_SYMBOLIC_VALUES);
  table_score_symb=NULL;

  free_multiregr_table_score();

  size_current_tree_table=0;
  size_current_tree_table_pred=0;

}

int allocate_tree_tables(int nb_of_nodes, int nb_of_leaves, int nb_pred, int tsp) {
  float *temp;

  free_tree_tables();

  left_successor=(int *)MyMalloc((size_t)nb_of_nodes*sizeof(int));
  if (left_successor==NULL) {
    return 0;
  }
  right_successor=(int *)MyMalloc((size_t)nb_of_nodes*sizeof(int));
  if (right_successor==NULL) {
    free_tree_tables();
    return 0; 
  }
  tested_attribute=(int *)MyMalloc((size_t)nb_of_nodes*sizeof(int));
  if (tested_attribute==NULL) {
    free_tree_tables();
    return 0;
  }

  node_weight=(float *)MyMalloc((size_t)nb_of_nodes*sizeof(float));
  if (node_weight==NULL) {
    free_tree_tables();
    return 0;
  }
  node_size=(float *)MyMalloc((size_t)nb_of_nodes*sizeof(float));
  if (node_size==NULL) {
    free_tree_tables();    
    return 0;
  }
  threshold=(union threshold_type *)MyMalloc((size_t)nb_of_nodes*sizeof(union threshold_type));
  if (threshold==NULL) {
    free_tree_tables();    
    return 0;
  }

  /* ne sert a rien si multiregr_savepred est a 1 */
  prediction=(int *)MyMalloc((size_t)nb_of_nodes*sizeof(int));
  if (prediction==NULL) {
    free_tree_tables();
    return 0;
  }

  nb_of_predictions=nb_pred;

  if (nb_pred>0) {

    /* printf("Table prediction_values = %d x %d\n",nb_of_leaves,nb_pred); */
    fflush(stdout);

    prediction_values=allocate_table_float(nb_of_leaves,nb_pred);
    if (prediction_values==NULL) {
      free_tree_tables();
      return 0;
    }
    
    /* allocation de la table de score (a ne pas faire si multiregr) */
    if (tsp==1) {
      table_score=allocate_table_score_type(3,nb_pred+1);
      if (table_score==NULL) {
	free_tree_tables();
	return 0;
      }
      table_score_symb=allocate_table_score_type(MAX_NUMBER_OF_SYMBOLIC_VALUES,nb_pred+1);
      if (table_score==NULL) {
	free_tree_tables();
	return 0;
      }
    }

    size_current_tree_table_pred=nb_of_leaves;

  } else
    size_current_tree_table_pred=0;

  size_current_tree_table=nb_of_nodes;

  return 1;

}

/* fonction qui recupere la taille du modele dernierement construit */

int get_last_model_complexity () {
  return (index_nodes+1);
}

/* set print_result */

void set_print_result(int i, int j) {
  print_result=i;
  print_detailled_result=j;
}

/* set best_first */

void set_best_first(int bf_flag, int tbf_flag, int max_tests) {
  best_first=bf_flag;
  true_best_first=tbf_flag;
  best_first_max_nb_tests=max_tests;

  if (size_tbf>0) {
    MyFree((int *)tbf_best_attribute);
    MyFree((int *)tbf_stop_splitting);
    MyFree((int *)tbf_not_significant);
    MyFree((SCORE_TYPE *)tbf_best_threshold_score);
    MyFree((SCORE_TYPE *)tbf_best_threshold_info);
    MyFree((union threshold_type  *)tbf_best_threshold);
    size_tbf=0;
  }
  
  if (tbf_flag) {
    tbf_best_attribute=(int *)MyMalloc((size_t)(max_tests+1)*sizeof(int));
    tbf_stop_splitting=(int *)MyMalloc((size_t)(max_tests+1)*sizeof(int));
    tbf_not_significant=(int *)MyMalloc((size_t)(max_tests+1)*sizeof(int));
    tbf_best_threshold_score=(SCORE_TYPE *)MyMalloc((size_t)(max_tests+1)*sizeof(SCORE_TYPE));
    tbf_best_threshold_info=(SCORE_TYPE *)MyMalloc((size_t)(max_tests+1)*sizeof(SCORE_TYPE));
    tbf_best_threshold=(union threshold_type *)MyMalloc((size_t)(max_tests+1)*sizeof(union threshold_type));    
    size_tbf=max_tests+1;
  }
}


/* fonction de test par defaut */

float *(*get_tree_prediction_vector)(int tree, int obj)=get_tree_prediction_vector_classical;

/* pour changer la fonction de test */

void set_test_classical() {
  get_tree_prediction_vector=get_tree_prediction_vector_classical;
}

/* recupere des elements des arbres. Permet de le reconstruire du côté lisp */

int get_tested_attribute(int node) {
  return tested_attribute[node];
}

int get_left_successor(int node) {
  if (left_successor[node]==-1)
    return -1;
  else
    return node+left_successor[node];
}

int get_right_successor(int node) {
  return node+right_successor[node];
}

float get_num_att_threshold(int node) {
  return threshold[node].f;
}

void get_symb_att_subset(int node, int *boolean_vector) {
  int i;
  for (i=0; i<NB_VAL(tested_attribute[node]); i++) {
    boolean_vector[i]=BITN(threshold[node],i);
  }
}

/*********************/
/* ENSEMBLE D'ARBRES */
/*********************/

/* description de l'ensemble d'arbres courant */
#define MAX_NUMBER_OF_TREES 10000
int ltrees[MAX_NUMBER_OF_TREES];
float ltrees_weight[MAX_NUMBER_OF_TREES];
int current_nb_of_ensemble_terms=0;
int average_predictions_ltrees=1; /* 1 -> divide pred by the sum of weights 0-> no */

/* Différentes variantes pour la creation du LS */

/* identite */

float make_ls_vector_identity(int tree) {
  /* on ne fait rien */
  return 1.0;
}

/* bagging */

float make_ls_vector_bagging(int tree) {
  int i,j;

  /* on remet les poids a zero */
  for (i=0; i<global_learning_set_size; i++)
    object_weight[i]=0.0;

  /* on incremente les poids d'objets tires au hasard */
  for (i=0; i<global_learning_set_size; i++) {
    int rn=get_random_integer(global_learning_set_size);
    object_weight[rn]+=1.0;
  }
  
  /* on construit le ls en prenant les objets de poids non nuls */
  current_learning_set_size=0;
  for (i=0; i<global_learning_set_size; i++) {
    if (object_weight[i]!=0.0) {
      current_learning_set[current_learning_set_size]=i;
      current_learning_set_size++;
    }
  }
  
  return 1.0;
}


/*----------------------------*/
/* MART regression (Friedman) */
/*----------------------------*/

/* différentes variantes pour le test des arbres */

float current_sum_weight=0.0;
float mart_mu=0.1;

void set_ensemble_method_parameters(int i, int nb_terms, int se, int sewg, float mu) {
  if (i==1) {
    make_ls_vector=make_ls_vector_bagging;
    average_predictions_ltrees=1;
  } else if (i==2) {
    make_ls_vector=make_ls_vector_mart_multiregr;
    mart_mu=mu;
    average_predictions_ltrees=0;
  } else {
    make_ls_vector=make_ls_vector_identity;
    average_predictions_ltrees=1;
  }

  store_ensemble=se;
  save_ensemble_while_growing=sewg;

  number_of_ensemble_terms=nb_terms;  
}

/* pour sauvegarder les ls */
int save_ensemble_ls=0;
int *save_ensemble_ls_vector=NULL;
float *save_ensemble_ls_weight=NULL;
int save_ensemble_ls_size[MAX_NUMBER_OF_TREES];
int save_ensemble_ls_pos=0;
int save_ensemble_ls_nb_ls=0;

void init_save_ensemble_ls(int b) {

  if (save_ensemble_ls && (save_ensemble_ls_vector!=NULL)) {
    MyFree((int *)save_ensemble_ls_vector);
    save_ensemble_ls_vector=NULL;
    MyFree((float *)save_ensemble_ls_weight);
    save_ensemble_ls_weight=NULL;
  }

  save_ensemble_ls=b;

  if (save_ensemble_ls) {
    int s=number_of_ensemble_terms*global_learning_set_size;
    save_ensemble_ls_vector=(int *)MyMalloc((size_t)s*sizeof(int));
    save_ensemble_ls_weight=(float *)MyMalloc((size_t)s*sizeof(float));
    save_ensemble_ls_pos=0;
    save_ensemble_ls_nb_ls=0;
  }
}

/* Construction de l'ensemble d'arbres */

float build_one_tree_ensemble(int *ts_vector, int length_ts_vector) {
  int i, j, t;
  int sum_complexity=0;
  float sum_error_single=0.0;
  float min_error; int index_min_error;
  int current_tree=-1; float current_weight; 
  float previous_weight;
  int nbn=0;
  FILE *fp;

  /* on vide tout */
  clean_all_trees();
  current_nb_of_ensemble_terms=0;
  current_sum_weight=0.0;

  /* verification de la memoire */
  if (best_first)
    nbn=2*best_first_max_nb_tests+1;
  else
    nbn=(2*global_learning_set_size-1);

  if (!store_ensemble) { /* ensemble non stocke */
    if (size_current_tree_table < nbn) {
      printf("memoire trop faible pour construire un arbre\n");
      return -1.0;
    }
  } else {/* ensemble stocke en memoire totalement */
    
    if (size_current_tree_table <(number_of_ensemble_terms*nbn)) {
      printf("memoire trop faible pour construire l'ensemble d'arbres\n");
      return -1.0;
    }
    if (number_of_ensemble_terms>MAX_NUMBER_OF_TREES) {
      printf("maximum %d trees, increase table size\n",MAX_NUMBER_OF_TREES);
      return -1.0;
    } 
  }

  if (save_ensemble_while_growing) {
    fp=fopen("temp-ensemble-trees.dat", "wb");
    fwrite(&average_predictions_ltrees, sizeof(float), 1, fp);
  }
  
  /* initialisation du test */
  if (print_result && length_ts_vector>0) {
    printf("initialisation de l'ensemble de test...\n");
    fflush(stdout);
  }
  /* allocation de la matrice de test si necessaire */
  /* de la taille (length_ts_vector*nb_pred) */
  
  /* initialisation de l'ensemble d'apprentissage */
  make_ls_vector(-1);
  
  /* boucle de construction */
  for (t=0; t<number_of_ensemble_terms; t++) {
    int current_tree;

    /* si on le demande, on sauve les LS */
    /* pour compute_node_subset_current_ensemble. Idealement, on devrait l'implementer
     * aussi pour le boosting en tenant compte des poids. Ca permettrait de faire des calculs
     * de variable importance plus precise.
     */
    if (save_ensemble_ls) {
      save_ensemble_ls_size[save_ensemble_ls_nb_ls]=current_learning_set_size;
      save_ensemble_ls_nb_ls++;
      for (i=0; i<current_learning_set_size; i++) {
	save_ensemble_ls_vector[save_ensemble_ls_pos]=current_learning_set[i];
	save_ensemble_ls_weight[save_ensemble_ls_pos]=object_weight[current_learning_set[i]];
	save_ensemble_ls_pos++;
      }
    }

    if (print_result) {
      printf("model %d...",t+1);
      fflush(stdout);
    }

    /* melange les attributes: */
    /* randomize_vector(attribute_vector, old_nb_attributes); */
    /* nb_attributes=nb_effective_att; */

    /* construction du modele */
    current_tree=build_one_tree();

    sum_complexity+=index_nodes-current_tree+1;
    current_weight=make_ls_vector(current_tree);

    if (print_result) {
      printf("(compl=%d) ", index_nodes-current_tree+1);
      fflush(stdout);
    }

    /* stockage du modele si on ne teste pas tout de suite */
    if (save_ensemble_while_growing) {
      fwrite(&current_weight,sizeof(float),1,fp);
      write_one_tree(current_tree, fp);
    }

    if (store_ensemble) {
      ltrees[t]=current_tree;
      ltrees_weight[t]=current_weight;
      current_nb_of_ensemble_terms++;
    } else
      clean_all_trees();

    if (print_result) printf("\n");
    
    if (current_weight==0.0) {
      /* on arrete, le modele precedent est trop mauvais */
      if (print_result) {
	printf("Stop (error too high)\n");
	fflush(stdout);
      } 
      t=number_of_ensemble_terms;
    } else if (current_weight<0) {
      /* on arrete, le dernier modele est parfait. Il devient le seul
	 modele */
      if (print_result) {
	printf("Stop (null error)\n");
	fflush(stdout);
      }
      t=number_of_ensemble_terms;
      if (store_ensemble) {
	current_nb_of_ensemble_terms=1;
	ltrees[0]=ltrees[t];
	ltrees_weight[0]=1.0;
      }
      sum_complexity=index_nodes-current_tree+1;
    }
  }
  if (print_result) {
    printf("complexity: %d",sum_complexity);
    printf("\n");
    fflush(stdout);
  }
  
  if (save_ensemble_while_growing) {
    fclose(fp);
    fp=fopen("temp-ensemble-nb-trees.dat", "wb");
    fwrite(&current_nb_of_ensemble_terms, sizeof(int), 1, fp);
    fclose(fp);
  }

  return -1.0;
}

/* pour pouvoir tester les premiers arbres d'un ensemble */

void set_current_nb_of_ensemble_terms(int nb) {
  current_nb_of_ensemble_terms=nb;
}

/* pour voir le nombre d'arbres construits */
int get_current_nb_of_ensemble_terms() {
  return current_nb_of_ensemble_terms;
}

/* pour pouvoir tester un arbre comme un ensemble d'un seul arbre */

void set_one_tree_ensemble() {
  current_nb_of_ensemble_terms=1;
  ltrees[0]=0;
  ltrees_weight[0]=1.0;
}

/* pour virer l'ensemble */

void clean_all_ensemble () {
  current_nb_of_ensemble_terms=0;
}

/* renvoit au c la racine de l'arbre i de l'ensemble */

int get_root_tree_ensemble(int i) {
  return ltrees[i];
}

float get_weight_tree_ensemble(int i) {
  return ltrees_weight[i];
}

void pretty_print_tree_from_ensemble(int tree, int max_depth) {
  prettyprinttree(ltrees[tree],max_depth);
}

int get_tree_leaf_from_ensemble(int tree, int obj) {
  return get_tree_leaf(ltrees[tree],obj);
}

int get_tree_pred_nb_from_ensemble(int tree, int obj) {
  return prediction[get_tree_leaf(ltrees[tree],obj)];
}

/****************************************
 * SAUVEGARDE ET CHARGEMENT DES MODÈLES *
 ****************************************/

/* calcule le nombre de noeud dans un arbre */

int get_tree_nb_nodes(int tree) {
  int nb_nodes=1;
  int current_node;
  
  index_stack_open_nodes=-1;
  index_stack_open_nodes++;
  stack_open_nodes[index_stack_open_nodes][0]=tree;
  
  while(index_stack_open_nodes>=0) {
    current_node=stack_open_nodes[index_stack_open_nodes][0];
    index_stack_open_nodes--;
    if (left_successor[current_node]!=-1) {
      nb_nodes+=2;
      index_stack_open_nodes++;
      stack_open_nodes[index_stack_open_nodes][0]=current_node+left_successor[current_node];
      index_stack_open_nodes++;
      stack_open_nodes[index_stack_open_nodes][0]=current_node+right_successor[current_node];
    }
  }

  return nb_nodes;

}

/* Ecrit un arbre dans le fichier (version propre a l'algorithme de construction) */

void write_one_tree(int tree, FILE *fp) {
  int nb_nodes=get_tree_nb_nodes(tree);
  int current_node;

  printf("nombre de noeuds a ecrire=%d\n", nb_nodes);
  fflush(stdout);
  /* on ecrit le nombre de noeuds */
  fwrite(&nb_nodes, sizeof(int), 1, fp);

  /* on ecrit l'info sur les noeuds */
  for (current_node=tree; current_node<tree+nb_nodes; current_node++) {
    int pred;
    /* on ecrit les info sur ce noeud dans le fichier */
    /* on ecrit 0 ou 1 selon la valeur de prediction[current_node] */
    /* deux cas: */
    if (left_successor[current_node]!=-1) {
      pred=-1;
      fwrite(&pred, sizeof(int), 1, fp);
      fwrite(&tested_attribute[current_node], sizeof(int), 1, fp);
      fwrite(&threshold[current_node], sizeof(union threshold_type), 1, fp);
      fwrite(&left_successor[current_node], sizeof(int), 1, fp);
      /* normalement, on n'a pas besoin de cette valeur */
      fwrite(&right_successor[current_node], sizeof(int), 1, fp);
    } else {
      pred=1;
      fwrite(&pred, sizeof(int), 1, fp);
      /* on ecrit les valeurs de predictions */
      fwrite(prediction_values[prediction[current_node]], sizeof(float), nb_classes, fp);
    }
  }
}

void read_one_tree(FILE *fp) {
  int nb_nodes;
  int i;
  
  /* on lit le nombre de noeuds */
  fread(&nb_nodes, sizeof(int), 1, fp);
  
  printf("nombre de noeuds a lire=%d\n", nb_nodes);
  fflush(stdout);
  
  if (index_nodes+nb_nodes>size_current_tree_table) {
    printf("not enough memory to load the tree (read_one_tree)\n");
    fflush(stdout);
  }
  /* on lit l'info sur tous les noeuds */
  for (i=0; i<nb_nodes; i++) {
    int pred;
    /* on rajoute un noeud */
    index_nodes++;
    fread(&pred, sizeof(int), 1, fp);
    if (pred==-1) {
      /* c'est un noeud test */
      prediction[index_nodes]=-1;
      fread(&tested_attribute[index_nodes], sizeof(int), 1, fp);
      fread(&threshold[index_nodes], sizeof(union threshold_type), 1, fp);
      fread(&left_successor[index_nodes], sizeof(int), 1, fp);
      fread(&right_successor[index_nodes], sizeof(int), 1, fp);
    } else {
      /* c'est une feuille */
      index_prediction_values++;
      prediction[index_nodes]=index_prediction_values;
      fread(prediction_values[index_prediction_values], sizeof(float), nb_classes, fp);
      tested_attribute[index_nodes]=-1;
      left_successor[index_nodes]=-1;
      right_successor[index_nodes]=-1;
    }
  }
}

void write_tree_to_file(int tree, char *filename) {
  FILE *fp=fopen(filename,"wb");
  write_one_tree(tree, fp);
  fclose(fp);
  printf("done");
  fflush(stdout);
}

void read_tree_from_file(char *filename) {
  FILE *fp=fopen(filename,"rb");
  clean_all_trees();
  read_one_tree(fp);
  fclose(fp);
  printf("done");
  fflush(stdout);
}

void write_ensemble_to_file(char *filename) {
  char str1[100];
  char str2[100];
  FILE *fp;
  FILE *fp_nb_trees;
  
  int i;

  strcpy(str1,filename);
  strcat(str1,"-trees.dat");
  fp=fopen(str1, "wb");
  strcpy(str2,filename);
  strcat(str2,"-nb-trees.dat");
  fp_nb_trees=fopen(str2, "wb");
  
  fwrite(&current_nb_of_ensemble_terms, sizeof(int), 1, fp_nb_trees);
  fclose(fp_nb_trees);

  /* on ecrit d'abord le nombre d'arbres, puis le poids des arbres */
  /* fwrite(&current_nb_of_ensemble_terms, sizeof(int), 1, fp); */
  fwrite(&average_predictions_ltrees, sizeof(int), 1, fp);
  /* fwrite(ltrees_weight, sizeof(float), current_nb_of_ensemble_terms, fp); */
  /* ensuite on ecrit les arbres un par un */
  for (i=0; i<current_nb_of_ensemble_terms; i++) {
    fwrite(&ltrees_weight[i],sizeof(float),1,fp);
    write_one_tree(ltrees[i], fp);
  }
  fclose(fp);
  
  printf("done");
  fflush(stdout);
}

void read_ensemble_from_file(char *filename) {
  char str1[100];
  char str2[100];
  FILE *fp;
  FILE *fp_nb_trees;
  
  int i;

  strcpy(str1,filename);
  strcat(str1,"-trees.dat");
  fp=fopen(str1, "rb");
  strcpy(str2,filename);
  strcat(str2,"-nb-trees.dat");
  fp_nb_trees=fopen(str2, "rb");

  clean_all_trees();
  
  /* on ecrit d'abord le nombre d'arbres, puis le poids des arbres */
  fread(&current_nb_of_ensemble_terms, sizeof(int), 1, fp_nb_trees);
  if (current_nb_of_ensemble_terms> MAX_NUMBER_OF_TREES) {
     printf("not enough memory to load the ensemble of trees (read_ensemble_from_file)\n");
     fflush(stdout);
     fclose(fp); fclose(fp_nb_trees);
     return;
  }

  fread(&average_predictions_ltrees, sizeof(int), 1, fp);
  /* fread(ltrees_weight, sizeof(float), current_nb_of_ensemble_terms, fp); */
  /* ensuite on ecrit les arbres un par un */
   printf("nb trees a lire=%d average_pred=%d\n", current_nb_of_ensemble_terms, average_predictions_ltrees);
  fflush(stdout);
  
  for (i=0; i<current_nb_of_ensemble_terms; i++) {
    fread(&ltrees_weight[i],sizeof(float), 1, fp);
    ltrees[i]=index_nodes+1;
    read_one_tree(fp);
  }

  fclose(fp);
  fclose(fp_nb_trees);
  
  printf("done");
  fflush(stdout);
}

void save_ensemble(char *filename) {
  write_ensemble_to_file(filename);
}

void load_ensemble(char *filename) {
  read_ensemble_from_file(filename);
}
