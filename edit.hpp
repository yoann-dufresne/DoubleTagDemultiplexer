

using namespace std;


/* Compute the edit distance matrix between s1 and s2
 * /!\ MAX word length = 250
 */
void compute_matrix(const string & s1, const string & s2);
/* Get the string representing the matrix for edit distance
 * /!\ The matrix must be computed before this function call
 */
string matrix_2_string ();
/* Get the edit distance
 * /!\ The matrix must be computed before this function call
 */
uint get_distance ();

/* Get a representation of the alignment only (without s1 and s2) */
string get_align_symbols (const string & s1, const string & s2);
/* Get the alignment of the strings */
string get_alignment (const string s1, const string s2);

