#ifndef H_SEQUENTIAL
#define H_SEQUENTIAL

#define SEQ1_MAX_LEN 3000
#define SEQ2_MAX_LEN 2000

void readInputFile(char *main_sequence, char ***sub_sequences, int *sub_sequences_count);
void findBestScore(char *main_seq, char *sub_seq, int score_matrix[26][26], int *score, int *offset, int *mutant);
void sequential(int argc, char **argv);
char mutant_char(char c);
void strToUpper(char *str);
char toUpper(char c);
#endif // H_SEQUENTIAL