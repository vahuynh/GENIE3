
/*
 * Author: Pierre Geurts (p.geurts@ulg.ac.be)
 * (c) 2002-2010
 *
 * (Please do not distribute without permission of the author)
 *
 */

/* MACRO, CONSTANTES */

#define VALUE_PI 3.141592653589793
#define PLUS_INFINI 1e100
#define MOINS_INFINI -1e100

/* precision des differentes fonction de calcul */

#define PRECISION double

#define MAX_OPEN_NODES 15000
#define MAX_TEST_CASES 20000

#define SCORE_TYPE double
/* permet de definir la core_table en char pour les grosses BDs */
#define CORETABLE_TYPE float
/* #define CORETABLE_TYPE unsigned char */

#define ERROR_TYPE double
#define MAX_NUMBER_OF_ATTRIBUTES 300000

/* pour les variables symboliques */

#define MAX_NUMBER_OF_SYMBOLIC_VALUES 1024
#define MAX_NUMBER_OF_SYMBOLIC_VALUES_DIV_32 32 /* C'est à dire maximum 128 valeurs */

#define BITN(x,n) (x.i[n/32]>>(n%32))%2
#define SET_BITN(x,n) x.i[n/32]|=(1<<(n%32))

#define SYMBOLICp(att) attribute_descriptors[att]>0
#define NUMERICALp(att) attribute_descriptors[att]==0
#define NB_VAL(att) attribute_descriptors[att]

#define MAX(x,y) (x<y? y: x)
#define MIN(x,y) (x<y? x: y)

/* redefinition du malloc (pour l'interface matlab) */
#define MyMalloc malloc
#define MyFree free

/* definition de types */

union threshold_type {
  unsigned int i[MAX_NUMBER_OF_SYMBOLIC_VALUES_DIV_32];
  float f;
};

/* prototype de fonctions */

/* fonction qui calcule des stats pour le vecteur et qui place le resultat dans table_score */
void (*summarize_vector)(int *vector, int start, int end);
/* fonction qui cree un vecteur de prediction */
int (*make_leaf_prediction)();
/* fonction qui renvoie true si on doit arreter la construction */
int (*stop_splitting_criterio)();
/* fonction d'evaluation pour le best_first */
SCORE_TYPE (*best_first_node_priority_function)();
/* fonction qui renvoie true si le test choisi n'est pas significatif */
int (*not_significant_test)();
/* fonction qui calcule le score a partir d'une table de score */
SCORE_TYPE (*compute_score_from_table)();
/* une fonction qui recherche un test pour splitter le noeud */
void (*find_a_split)(int *ls_vector, int start, int end);
/* une fonction qui recherche un seuil pour un attribut numerique */
void (*find_a_threshold_num)(int att, int *ls_vector, int start, int end);
/* une fonction qui recherche un seuil pour un attribut symbolique */
void (*find_a_threshold_symb)(int att, int *ls_vector, int start, int end);

void write_one_tree(int tree, FILE *fp);

/* Sorting function */

void quicksort_ls_vector(int ls_vector[], int start, int end, int att);
int separate_ls_vector_local(int best_attribute, union threshold_type best_threshold, int *ls_vector, int start, int end);

/* generique */
int build_one_tree();
void find_a_split_at_random(int *ls_vector, int start, int end);
void find_a_threshold(int att, int *ls_vector, int start, int end);
int check_test(int att, float val, union threshold_type threshold);

/* utilitaires */
int get_random_integer(int max_val);
float get_random_normal();
void get_random_subset_non_empty_full(union threshold_type *rs, int max_val);
float get_random_float();
int my_best_class(float *proba_vector, int nb_classes);
float chi_square(int ddl, float conf);

float getattval(int obj, int att);

/* interface */
void clean_all_trees();
void clean_all_ensemble();
void set_random_threshold_method(int i);

float make_ls_vector_bagging(int tree);
float make_ls_vector_identity(int tree);

int get_tree_nb_nodes(int tree);

void free_multiregr_table_score();
float make_ls_vector_mart_multiregr(int tree);
