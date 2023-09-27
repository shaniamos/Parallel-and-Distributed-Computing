#ifndef H_SEQUENTIAL
#define H_SEQUENTIAL

#define SEQ1_MAX_LEN 3000
#define SEQ2_MAX_LEN 2000
enum
{
    RESULT_SCORE,
    RESULT_OFFSET,
    RESULT_MUTANT
};

void readInputFile(char *main_sequence, char ***sub_sequences, int *sub_sequences_count);
void findBestScore(char *main_seq, char *sub_seq, int score_matrix[26][26], int *score, int *offset, int *mutant, int rank, int num_procs);
void parallel(int argc, char *argv[], int rank, int num_procs);
char mutate_char(char c);
void strToUpper(char *str);
char toUpper(char c);
void cuda_init(char *main_seq, char *sub_seq, int score_matrix[26][26], int main_seq_len, int sub_seq_len, char **dev_main_seq, char **dev_sub_seq, int **dev_score_matrix, int **dev_best_score);

void cuda_free(char *dev_main_seq, char *dev_sub_seq, int *dev_score_matrix, int *dev_best_scores);
void findBestOffsetCuda(char *dev_main_seq, char *dev_sub_seq, int *dev_score_matrix, int *dev_best_scores, int *best_score, int *best_offset, int *cuda_results, int mutant, int main_seq_len, int sub_seq_len);
#endif // H_SEQUENTIAL