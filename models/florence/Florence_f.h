/** florence_f.h **/

void florencef77_(double PINP[28], double *FREQ, double R1[10]);
void modflor_(double PINP[28], double *FREQ, double R1[10]);

/* C side, not a model function: RUNO(1..5) from one evaluation (Florence_c.c),
 * parameters in the PINP order of FlorenceN without index and n. */
void FlorenceRUNO(double FREQ, const double parameters[28], double runo[5]);
